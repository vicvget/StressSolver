#include "StressStrainFortranSolver.h"

#include "CsrSymmetricMatrix.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"

#include <omp.h>


using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;


//#define NO_LINKSH
//#define NO_INTOMSUB

#ifndef M_PI

#define M_PI 3.1415926535897932384626433832795

#endif

#define MAX(x, y) ((x) > (y) ? x : y)


double DMOD
	(
		const double d1,
		const double d2
	)
{
	return d1 - d2 * ((int)(d1 / d2));
}

StressStrainFortranSolver::StressStrainFortranSolver
	(
		double* params, 
		int* links, 
		int nLinks, 
		double *gridElements, 
		int nElements, 
		double gridStep, 
		double timeStep,
		int numThreads,
		int stride
	)
	:
		StressStrainSolver(nElements,stride),
		_isFirstSolution(true),
		_poissonRatio(0.)
{
	_isFirstIteration = true;
	_nIteration = 0;
	_nVariables = nElements * vecStride2;
	METS = 0;

	T = 0;										// время
	H = timeStep;								// шаг по времени
	_gridStep = gridStep;						// шаг сетки
	_elasticModulus = params[0];				// модуль упругости
	_dampingFactor = params[1];					// коэффициент демпфирования
	_density = params[2];						// плотность материала
	_stiffScale = params[3];					// масштабный коэффициент
	
	_cellMass = _density * _gridStep * _gridStep * _gridStep; // масса ячейки
	_elasticModulusScaled = _elasticModulus / _stiffScale;	  // отмасштабированный модуль упругости
	NumThreads = (numThreads > 0) ? numThreads : omp_get_max_threads();
	
	const int outWidth = 15;
	std::cout << "------------------------------" << std::endl
		<< "    FOR SOLVER IS CREATED" << std::endl
		<< "------------------------------" << std::endl
		<< std::setw(outWidth) << "VECSTRIDE: " << std::setw(outWidth) << vecStride << std::endl
		<< std::setw(outWidth) << "STIFF SCALE: " << std::setw(outWidth) << _stiffScale << std::endl
		<< std::setw(outWidth) << "ELASTIC: " << std::setw(outWidth) << _elasticModulus << std::endl
		<< std::setw(outWidth) << "DAMPING: " << std::setw(outWidth) << _dampingFactor << std::endl
		<< std::setw(outWidth) << "DENSITY: " << std::setw(outWidth) << _density << std::endl
		<< std::setw(outWidth) << "ELEMENTS: " << std::setw(outWidth) << nElements << std::endl
		<< std::setw(outWidth) << "GRID STEP: " << std::setw(outWidth) << _gridStep << std::endl;

	//GRA = _dataRotationMtx;//new double [nElements * 9];

	//for (int j = 0; j < nElements; j++)
	//{
	//	for (int i = 0; i < 9; i++)
	//	{
	//		_dataRotationMtx[j * 9 + i] = 0;
	//	}
	//	_dataRotationMtx[j * 9] = _dataRotationMtx[j * 9 + 4] = _dataRotationMtx[j * 9 + 8] = 1;
	//}


	_initX = new double[_nVariables];
	_initDX= new double[_nVariables];
	_hDDX1 = new double[_nVariables];
	_hDDX2 = new double[_nVariables];
	_hDDX3 = new double[_nVariables];

	_varX = _dataInternal;
	_varDX = _dataInternal+_nVariables;
	_varDDX = _dataInternal+_nVariables*2;

	R = new double[_nVariables];
	RZ = new double[_nVariables];
	R1Z = new double[_nVariables];
	_fue = new Fue(nElements);
	_copym = new Copym(nElements);
	MLINK = new int[nElements * 6];
	_nodes = new double[nElements*3];
	memcpy(_nodes, gridElements, sizeof(double) * nElements * 3);
	memset(MLINK, 0, nElements * 6 * sizeof(int));
	for (int i = 0; i < nLinks * 2; i += 2)
	{
		double* node1 = &_nodes[3 * links[i]];
		double* node2 = &_nodes[3 * links[i + 1]];
		for (int k = 0; k < 3; k++)
		{
			double dx = node1[k] - node2[k];
			if ((dx < 0) && (fabs(gridStep + dx) < gridStep * 1e-4))
			{
				MLINK[links[i + 1] * 6 + k] = links[i] + 1;
				MLINK[links[i] * 6 + k + 3] = links[i + 1] + 1;
			}
			else if (!(dx < 0) && (fabs(gridStep - dx) < gridStep * 1e-4))
			{
				MLINK[links[i + 1] * 6 + k + 3] = links[i] + 1;
				MLINK[links[i] * 6 + k] = links[i + 1] + 1;
			}
		}
	}
	memset(_dataInternal, 0, _nVariables * 3 * sizeof(double));

	//// единичные матрицы поворота
	//for (size_t i = 0; i < _nElements; i++)
	//{
	//	_dataRotationMtx[i * 9] = 1.;
	//	_dataRotationMtx[i * 9 + 3 + 1] = 1.;
	//	_dataRotationMtx[i * 9 + 6 + 2] = 1.;
	//}

	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_dataInternal[i * vecStride2 + j] = gridElements[i * 3 + j];
		}
	}
	intomsub();
	_testTimer.Allocate(10);
}

// virtual
StressStrainFortranSolver::~StressStrainFortranSolver()
{
	//delete [] GR1;
	//delete [] GRA;

	delete [] _initX;
	delete [] _initDX;
	delete [] _hDDX1;
	delete [] _hDDX2;
	delete [] _hDDX3;

	delete [] R;
	delete [] RZ;
	delete [] R1Z;
	delete [] MLINK;
	delete [] _nodes;
	delete _fue;
	delete _copym;
}

/** Добавление связей
* @param links - массив 6*_nElements с признаками связей (1 - есть связь, 0 - нет связи)
*/
// virtual
void StressStrainFortranSolver::AddLinks
	(
		const int* links
	)
{
	for (int i = 0; i < _nElements * 6; i++)
	{
		MLINK[i] = links[i];
	}
}

void StressStrainFortranSolver::AddPartialBoundary
	(
		int* boundaryNodesIndicesInPart, 
		int numberOfNodesInPart,
		int numberOfNodes,
		int bcKind,
		double* bcParams
	)
{
	AddBoundary(boundaryNodesIndicesInPart,numberOfNodesInPart,bcKind,bcParams);
	_boundaryParamsSet.back().SetNodesCountInFullBoundary(numberOfNodes);
}

// virtual
void StressStrainFortranSolver::AddBoundary
	(
		int* boundaryNodesIndices, 
		int numberOfBoundaryNodes,
		int bcKind,
		double* bcParams
	)
{
	double bcParamsCopy[7];

	switch (bcKind)
	{
	case 3:
		for (int i = 0; i < 6; i++)
		{
			bcParamsCopy[i] = bcParams[i] / _stiffScale;
		}
		bcParamsCopy[6] = bcParams[6];
		break;
	
	case 4:
		for (int i = 0; i < 6; i++)
		{
			bcParamsCopy[i] = bcParams[i];
		}
		break;

	default:
		// TODO: вставить exception
		;	// Ошибка!!!
	}

	BoundaryParams bp
		(
			bcKind,
			bcParamsCopy,
			boundaryNodesIndices,
			numberOfBoundaryNodes
		);

	_boundaryParamsSet.push_back(bp);
}

// virtual
void StressStrainFortranSolver::ChangeBoundaryParam
	(
		const int bcNumber,
		const int bcParamNumber,
		const double bcParamValue
	)
{
	_boundaryParamsSet[bcNumber].SetParam
		(
			bcParamValue / _stiffScale,
			bcParamNumber
		);
}

// virtual
double StressStrainFortranSolver::GetBoundaryParam
	(
		const int bcNumber,
		const int bcParamNumber
	)
{
	return _boundaryParamsSet[bcNumber].GetParam(bcParamNumber);
}

// virtual
void StressStrainFortranSolver::UpdateBuffer
	(
		double scaleFactor
	)
{
	GetScalarParameter(_data);

	float* coordinatesData = _data + _nElements;

	GetVectorParameter(coordinatesData);
}

void StressStrainFortranSolver::intomsub()
{
	//return;
#ifdef NOINTOMSUB
	return;
#endif

	double* GR1 = _dataInternal;
	int i;

//#pragma omp parallel for private(i) num_threads(NumThreads)
	//std::cout << "                _nElements = " << _nElements << std::endl;
	for (i = 0; i < _nElements; i++)
	{
		urejl4s
			(
				&_dataRotationMtx[i * 9],
				&_dataInternal[i * 6 + 3],
				&_dataInternal[i * 6 + 3 + _nElements * 6],
				METS,
				i
			);
	}
}

void StressStrainFortranSolver::urejl4s
	(
		double* a,
		double *ug,
		double *om,
		int mets,
		const int id
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	int K = id * 3;
	int KA = id * 9;
	
	if (mets == 0) 
	{
		_fue->Update(K);
		_fue->FLAGRY[id] = 0;
		_copym->Copy(a, _fue->A1 + KA, id, T);
		
		//std::cout << "#############################################" << std::endl;
		//std::cout << "intomsub SolveElementRotation METS == 0" << std::endl;

		return;
	}
	
	double UY = _fue->UE[K];
	UY = fabs(DMOD(UY, 2 * M_PI));           // DMOD - ???????

	if ((mets == 1) && (UY > 1.28))
	{

		_fue->Update(K);
		_fue->FLAGRY[id] = 1.0;
		_fue->IFLAGRY = 1;
		
		_copym->Copy(a, _fue->A1 + KA, id, T);
	} 
	_fue->Update(K, mets);
	_fue->UpdateR(K, om, H);
	_fue->UpdateR2(K, mets);

	_fue->UpdateMtx(K, a);
	multm2(_fue->A1 + KA, a);
}

void StressStrainFortranSolver::multm2
	(
		double *a1,
		double *a2
	)
{
	double az[9];

	for (int i = 0; i < 9; i++)
	{
		az[i] = a2[i];
	}
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)  
		{
			// TODO: переписать через цикл
			a2[j * 3 + i] = a1[j * 3] * az[i] + a1[j * 3 + 1] * az[3 + i] + a1[j * 3 + 2] * az[6 + i];
		}
	}
}

void StressStrainFortranSolver::linksh
	(
		double *X1,
		double *X2,
		double *SL,
		double *VL,
		double *A,
		double* C,
		int j,
		int N2,
		int KUZ
	)
{	
	//std::cout << "!";
	double* GR1 = _dataInternal;
	double* GRA = _dataRotationMtx;

	double *A1 = &GRA[j * 9];
	double *A2 = &GRA[N2 * 9];

	double *P1 = &GR1[j * 6];
	double *P2 = &GR1[N2 * 6];

	double *V1 = &GR1[j * 6 + KUZ * 6];
	double *V2 = &GR1[N2 * 6 + KUZ * 6];

	double *OM1 = &GR1[j * 6 + KUZ * 6 + 3];
	double *OM2 = &GR1[N2 * 6 + KUZ * 6 + 3];


	double B[12];

	double EM0010 = A1[0];
	double EM0011 = A1[1];
	double EM0012 = A1[2];
	double EM0013 = A1[3];
	double EM0014 = A1[4];
	double EM0015 = A1[5];
	double EM0016 = A1[6];
	double EM0017 = A1[7];
	double EM0018 = A1[8];
	double EM0024 = OM1[0];
	double EM0025 = OM1[1];
	double EM0020 = OM1[2];
	double EM0028 = P1[0];
	double EM0029 = P1[1];
	double EM0030 = P1[2];
	double EM0031 = V1[0];
	double EM0032 = V1[1];
	double EM0033 = V1[2];
	double EM0034 = OM1[0];
	double EM0035 = OM1[1];
	double EM0036 = OM1[2];
	double EM0037 = X1[0];
	double EM0038 = X1[1];
	double EM0039 = X1[2];
	double EM0040 = A2[0];
	double EM0041 = A2[1];
	double EM0042 = A2[2];
	double EM0043 = A2[3];
	double EM0044 = A2[4];
	double EM0045 = A2[5];
	double EM0046 = A2[6];
	double EM0047 = A2[7];
	double EM0048 = A2[8];
	double EM0054 = OM2[0];
	double EM0055 = OM2[1];
	double EM0050 = OM2[2];
	double EM0058 = P2[0];
	double EM0059 = P2[1];
	double EM0060 = P2[2];
	double EM0061 = V2[0];
	double EM0062 = V2[1];
	double EM0063 = V2[2];
	double EM0064 = OM2[0];
	double EM0065 = OM2[1];
	double EM0066 = OM2[2];
	double EM0067 = X2[0];
	double EM0068 = X2[1];
	double EM0069 = X2[2];
	double EM0085 = EM0010;
	double EM0088 = EM0011;
	double EM0091 = EM0012;
	double EM0086 = EM0013;
	double EM0089 = EM0014;
	double EM0092 = EM0015;
	double EM0087 = EM0016;
	double EM0090 = EM0017;
	double EM0093 = EM0018;
	double EM0104 = EM0039;
	double EM0132 = EM0037;
	double EM0136 = EM0038;
	double EM0142 = P1[3];
	double EM0143 = P2[3];
	double EM0144 = P1[4];
	double EM0145 = P2[4];
	double EM0146 = P1[5];
	double EM0147 = P2[5];
	//  GENA
	//    произведение матрицы поворота второго тела на кооординаты узла
	double EM0070 = +EM0040 * EM0067 + EM0043 * EM0068 + EM0046 * EM0069;
	double EM0071 = +EM0041 * EM0067 + EM0044 * EM0068 + EM0047 * EM0069;
	double EM0072 = +EM0042 * EM0067 + EM0045 * EM0068 + EM0048 * EM0069;
	// VEKMIN
	//     разница координат ц.м. второго и первого тел
	double EM0073 =										+EM0058 - EM0028;
	double EM0074 =										+EM0059 - EM0029;
	double EM0075 =										+EM0060 - EM0030;
	// VEKMIN END
	// VEKSUM
	//     сумма предыдущих перемещений
	double EM0076 =										+EM0073 + EM0070;
	double EM0077 =										+EM0074 + EM0071;
	double EM0078 =										+EM0075 + EM0072;
	// VEKSUM END
	//     умножение транспонированной матрицы первого тела на предыдущее перемещение
	double EM0079 = +EM0010 * EM0076 + EM0011 * EM0077 + EM0012 * EM0078;
	double EM0080 = +EM0013 * EM0076 + EM0014 * EM0077 + EM0015 * EM0078;
	double EM0081 = +EM0016 * EM0076 + EM0017 * EM0077 + EM0018 * EM0078;
	// VEKMIN
	//     вычитание из координат ц.м. первого тела координат узла
	double EM0082 =                            +EM0079-EM0037;
	double EM0083 =                            +EM0080-EM0038;
	double EM0084 =                            +EM0081-EM0039;
	// VEKMIN END
	//     произведение элементов матриц поворота
	double EM0094 =+EM0085*EM0040+EM0088*EM0041+EM0091*EM0042;
	double EM0097 =+EM0085*EM0043+EM0088*EM0044+EM0091*EM0045;
	double EM0100 =+EM0085*EM0046+EM0088*EM0047+EM0091*EM0048;
	double EM0095 =+EM0086*EM0040+EM0089*EM0041+EM0092*EM0042;
	double EM0098 =+EM0086*EM0043+EM0089*EM0044+EM0092*EM0045;
	double EM0101 =+EM0086*EM0046+EM0089*EM0047+EM0092*EM0048;
	double EM0096 =+EM0087*EM0040+EM0090*EM0041+EM0093*EM0042;
	double EM0099 =+EM0087*EM0043+EM0090*EM0044+EM0093*EM0045;
	double EM0102 =+EM0087*EM0046+EM0090*EM0047+EM0093*EM0048;
	A[0]=EM0085;
	A[1]=EM0088;
	A[2]=EM0091;

	VL[0]=V1[0]*EM0085;
	VL[0]+=V1[1]*EM0088;
	VL[0]+=V1[2]*EM0091;
	double EM0105=                                   -EM0038;
	VL[0]+=OM1[1]*EM0104;
	VL[0]+=OM1[2]*EM0105;
	VL[0]-=V2[0]*EM0085;
	VL[0]-=V2[1]*EM0088;
	VL[0]-=V2[2]*EM0091;
	
	C[0]=-EM0085;
    C[1]=-EM0088;
    C[2]=-EM0091;
	
	A[3]=0.0;
	A[4]=EM0104;
	A[5]=EM0105;
	double EM0106=              +EM0100*EM0068-EM0097*EM0069;
	double EM0107=              +EM0094*EM0069-EM0100*EM0067;
	double EM0108=              +EM0097*EM0067-EM0094*EM0068;
	VL[0]-=OM2[0]*EM0106;
	VL[0]-=OM2[1]*EM0107;
	VL[0]-=OM2[2]*EM0108;
	C[3]=-EM0106;
	C[4]=-EM0107;
	C[5]=-EM0108;
	double EM0109=              -EM0020*EM0038+EM0025*EM0039;
	double EM0110=              +EM0020*EM0037-EM0024*EM0039;
	double EM0111=              -EM0025*EM0037+EM0024*EM0038;
	double EM0112=              -EM0020*EM0110+EM0025*EM0111;
	double EM0113=              +EM0020*EM0109-EM0024*EM0111;
	double EM0114=              -EM0025*EM0109+EM0024*EM0110;
	//**********************************
	double EM0115=EM0112;
	double EM0116=EM0113;
	double EM0117=EM0114;
	//**********************************
	double EM0118=              -EM0050*EM0068+EM0055*EM0069;
	double EM0119=              +EM0050*EM0067-EM0054*EM0069;
	double EM0120=              -EM0055*EM0067+EM0054*EM0068;
	double EM0121=              -EM0050*EM0119+EM0055*EM0120;
	double EM0122=              +EM0050*EM0118-EM0054*EM0120;
	double EM0123=              -EM0055*EM0118+EM0054*EM0119;
	double EM0124=+EM0094*EM0121+EM0097*EM0122+EM0100*EM0123;
	double EM0125=+EM0095*EM0121+EM0098*EM0122+EM0101*EM0123;
	double EM0126=+EM0096*EM0121+EM0099*EM0122+EM0102*EM0123;
	// VEKMIN
	double EM0127=                            +EM0124-EM0115;
	double EM0128=                            +EM0125-EM0116;
	double EM0129=                            +EM0126-EM0117;
	// VEKMIN END
	B[0] = +EM0127;
	SL[0] = -EM0082;

	VL[1]=V1[0]*EM0086;
	VL[1]+=V1[1]*EM0089;
	VL[1]+=V1[2]*EM0092;
	A[6]=EM0086;
	A[7]=EM0089;
	A[8]=EM0092;
	double EM0130=                                   -EM0039;
	VL[1]+=OM1[0]*EM0130;
	VL[1]+=OM1[2]*EM0132;
	VL[1]-=V2[0]*EM0086;
	VL[1]-=V2[1]*EM0089;
	VL[1]-=V2[2]*EM0092;
	A[9]=EM0130;
	A[10]=0.0;
	A[11]=EM0132;

	C[6]=-EM0086;
	C[7]=-EM0089;
	C[8]=-EM0092;
	double EM0133=              +EM0101*EM0068-EM0098*EM0069;
	double EM0134=              +EM0095*EM0069-EM0101*EM0067;
	double EM0135=              +EM0098*EM0067-EM0095*EM0068;
	VL[1]-=OM2[0]*EM0133;
	VL[1]-=OM2[1]*EM0134;
	VL[1]-=OM2[2]*EM0135;
	B[1] = +EM0128;
	SL[1] = -EM0083;

	C[9]=-EM0133;
	C[10]=-EM0134;
	C[11]=-EM0135;

	VL[2]=V1[0]*EM0087;
	VL[2]+=V1[1]*EM0090;
	VL[2]+=V1[2]*EM0093;
	A[12]=EM0087;
	A[13]=EM0090;
	A[14]=EM0093;
	double EM0137=                                   -EM0037;
	VL[2]+=OM1[0]*EM0136;
	VL[2]+=OM1[1]*EM0137;
	VL[2]-=V2[0]*EM0087;
	VL[2]-=V2[1]*EM0090;
	VL[2]-=V2[2]*EM0093;
	A[15]=EM0136;
	A[16]=EM0137;
	A[17]=0.0;
	C[12]=-EM0087;
	C[13]=-EM0090;
	C[14]=-EM0093;
	double EM0139=              +EM0102*EM0068-EM0099*EM0069;
	double EM0140=              +EM0096*EM0069-EM0102*EM0067;
	VL[2]-=OM2[0]*EM0139;
	VL[2]-=OM2[1]*EM0140;

	double EM0141=              +EM0099*EM0067-EM0096*EM0068;
	C[15]=-EM0139;
	C[16]=-EM0140;
	C[17]=-EM0141;

	VL[2]-=OM2[2]*EM0141;//*****************************************************
	B[2] = +EM0129;
	SL[2] = -EM0084;
	VL[3]=OM1[0]-OM2[0];

	for (int i=3;i<6;i++)
		for (int j=0;j<6;j++)
		{
			A[i*6+j]=0.0;
			C[i*6+j]=0.0;
		}

	A[21]=1.0;
	C[21]=-1.0;
	B[3]=0.0;
	SL[3]=EM0142-EM0143;
	VL[4]=OM1[1]-OM2[1];
	A[28]=1.0;
	C[28]=-1.0;
	B[4]=0.0;
	SL[4]=EM0144-EM0145;
	VL[5]=OM1[2]-OM2[2];
	A[35]=1.0;
	C[35]=-1.0;
	B[5]=0.0;
	SL[5]=EM0146-EM0147;

	//CompareAngles(A1, A2, SL);
}



void StressStrainFortranSolver::linksh3
	(
		double *cVec1,	//координаты точки связи элемента 1
		double *cVec2,	//координаты точки связи элемента 2
		double *SL,		// выход деформаций
		double *VL,		// выход изм. скоростей
		double *A,		// выход
		double* C,		// выход 
		int nodeId1,	// номер узла 1
		int nodeId2,	// номер узла 2
		int nNodes		// количество узлов
	)
{	
	const unsigned int strideVec = 3;							// смещение между векторами (линейные X,Y,Z или угловые Rx,Ry,Rz степени свободы)
	const unsigned int strideMat = strideVec * strideVec;		// смещение между матрицами поворота
	const unsigned int strideDeriv = nNodes * 2 * strideVec;	// смещение до производных
	const unsigned int nodeOffset1 = nodeId1 * 2 * strideVec;	// смещение до вектора неизвестных тела nodeId1
	const unsigned int nodeOffset2 = nodeId2 * 2 * strideVec;	// смещение до вектора неизвестных тела nodeId2

	Mat3 matA01(_dataRotationMtx + nodeId1 * strideMat); // смещение до матрицы поворота тела nodeId1
	Mat3 matA02(_dataRotationMtx + nodeId2 * strideMat); // смещение до матрицы поворота тела nodeId2
	Vec3Ref vecC1 = MakeVec3(cVec1);					 // координаты точки связи узла 1
	Vec3Ref vecC2 = MakeVec3(cVec2);					 // координаты точки связи узла 2
	
	// т.к. мы храним матрицу по строкам, а для fortran - хранится по столбцам
	matA02 = matA02.Tr();
	matA01 = matA01.Tr();
	Mat3 matA12 = matA01.Tmul(matA02);					 // A_{12} = A^T_{01}xA_{02}
	double *pVec1 = _dataInternal + nodeOffset1;	     // адрес смещение до вектора неизвестных тела nodeId1
	double *pVec2 = _dataInternal + nodeOffset2;		 // адрес смещение до вектора неизвестных тела nodeId2
	Vec3Ref vecP1 = MakeVec3(pVec1);					 // поступательные координаты тела nodeId1
	Vec3Ref vecP2 = MakeVec3(pVec2);					 // поступательные координаты тела nodeId2
	Vec3Ref vecR1 = MakeVec3(pVec1 + strideVec);		 // вращательные координаты тела nodeId1
	Vec3Ref vecR2 = MakeVec3(pVec2 + strideVec);		 // вращательные координаты тела nodeId2
	Vec3Ref vecV1 = MakeVec3(pVec1 + strideDeriv);       // линейные скорости тела nodeId1
	Vec3Ref vecV2 = MakeVec3(pVec2 + strideDeriv);		 // линейные скорости тела nodeId2
	Vec3Ref vecW1 = MakeVec3(pVec1 + strideDeriv + strideVec);   // угловые скорости тела nodeId1
	Vec3Ref vecW2 = MakeVec3(pVec2 + strideDeriv + strideVec);	 // угловые скорости тела nodeId2

	// переводим вектор линии точек связи С2-С1 в СК1
	Vec3 vecT0 = 
		vecC1 
		- matA12*vecC2 
		- matA01.Tmul(vecP2 - vecP1);
	
	vecT0.Export(SL);

	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) 
		+ vecW1.Cross(vecC1) 
		- matA12*(vecW2.Cross(vecC2));

	VecT1.Export(VL);

	// для угловых степеней свободы - просто берутся разница углов и угловых скоростей
	(vecR1 - vecR2).Export(SL + strideVec);
	(vecW1 - vecW2).Export(VL + strideVec);

	// Заполнение матрицы A
	//      a00 a01 a02 a03 a04 a05 
	//      a06 a07 a08 a09 a10 a11
	//      a12 a13 a14 a15 a16 a17
	//      a18 a19 a20 a21 a22 a23
	//      a24 a25 a26 a27 a28 a29
	//      a30 a31 a32 a33 a34 a35 
	//A
	//       A0  A1  A2   0  -X2 X1 
	//       A3  A4  A5  X2   0 -X0
	//       A6  A7  A8 -X1  X0   0
	//        0   0   0   1   0   0
	//        0   0   0   0   1   0
	//        0   0   0   0   0   1 

	for (int i=3;i<6;i++)
		for (int nodeId1=0;nodeId1<6;nodeId1++)
		{
			A[i*6+nodeId1]=0.0;
		}

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			A[i*6+j] = matA01[i*3+j];
		}
	}
	
	A[3]=0.0;        A[4]=-vecC1[2]; A[5]=vecC1[1]; 
	A[9]=vecC1[2];   A[10]=0.0;      A[11]=-vecC1[0];
	A[15]=-vecC1[1]; A[16]=vecC1[0]; A[17]=0.0;

	A[21]=1.0;
	A[28]=1.0;
	A[35]=1.0;
}

void CrossProduct
(
double* v1,
double* v2,
double* res
)
{
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = -v1[0] * v2[2] + v1[2] * v2[0];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void StressStrainFortranSolver::linksh2
	(
		double *cVec1,	//координаты точки связи элемента 1
		double *cVec2,	//координаты точки связи элемента 2
		double *SL,		// выход деформаций
		double *VL,		// выход изм. скоростей
		double *A,		// выход
		double* C,		// выход 
		int nodeId1,	// номер узла 1
		int nodeId2,	// номер узла 2
		int nNodes		// количество узлов
	)
{	
	if(nodeId1 == 100 && _nIteration == 10)
		int debug = 0;
	//std::cout << "!";
	//double* GR1 = _dataInternal;
	int strideSize = nNodes * 6;
	int nodeOffset1 = nodeId1 * 6;
	int nodeOffset2 = nodeId2 * 6;
	double* GRA = _dataRotationMtx;
	double *aMtx1 = GRA + nodeId1 * 9;
	double *aMtx2 = GRA + nodeId2 * 9;

	double *pVec1 = _dataInternal + nodeOffset1;
	double *pVec2 = _dataInternal + nodeOffset2;

	double *vVec1 = pVec1 + strideSize;
	double *vVec2 = pVec2 + strideSize;

	double *wVec1 = vVec1+3;
	double *wVec2 = vVec2+3;

	double *rVec1 = pVec1+3;
	double *rVec2 = pVec2+3;

	//  GENA
	//    произведение матрицы поворота второго тела на кооординаты узла A_{02}*cVec2
	double tVec1[3];
	tVec1[0] = +aMtx2[0] * cVec2[0] + aMtx2[3] * cVec2[1] + aMtx2[6] * cVec2[2];
	tVec1[1] = +aMtx2[1] * cVec2[0] + aMtx2[4] * cVec2[1] + aMtx2[7] * cVec2[2];
	tVec1[2] = +aMtx2[2] * cVec2[0] + aMtx2[5] * cVec2[1] + aMtx2[8] * cVec2[2];
	// VEKMIN
	//     разница координат ц.м. второго и первого тел pVec2- pVec1
	double tVec2[3];
	tVec2[0] =	pVec2[0] - pVec1[0];
	tVec2[1] =	pVec2[1] - pVec1[1];
	tVec2[2] =	pVec2[2] - pVec1[2];
	// VEKMIN END
	// VEKSUM
	//     сумма предыдущих перемещений
	double tVec3[3];
	tVec3[0] = tVec2[0] + tVec1[0];
	tVec3[1] = tVec2[1] + tVec1[1];
	tVec3[2] = tVec2[2] + tVec1[2];
	// VEKSUM END
	//     умножение транспонированной матрицы первого тела на предыдущее перемещение
	double tVec4[3];
	tVec4[0] = aMtx1[0] * tVec3[0] + aMtx1[1] * tVec3[1] + aMtx1[2] * tVec3[2];
	tVec4[1] = aMtx1[3] * tVec3[0] + aMtx1[4] * tVec3[1] + aMtx1[5] * tVec3[2];
	tVec4[2] = aMtx1[6] * tVec3[0] + aMtx1[7] * tVec3[1] + aMtx1[8] * tVec3[2];
	// VEKMIN
	//     вычитание из координат ц.м. первого тела координат узла
	double tVec5[3];
	tVec5[0] = tVec4[0]-cVec1[0];
	tVec5[1] = tVec4[1]-cVec1[1];
	tVec5[2] = tVec4[2]-cVec1[2];
	// VEKMIN END

	SL[0]=-tVec5[0];
	SL[1]=-tVec5[1];
	SL[2]=-tVec5[2];
	SL[3]=rVec1[0]-rVec2[0];
	SL[4]=rVec1[1]-rVec2[1];
	SL[5]=rVec1[2]-rVec2[2];

	//     произведение элементов матриц поворота
	// Зранение по столбцам A_{12}=A_{01}^T*A_{02}
	double tMtx[9];
	tMtx[0] = aMtx1[0]*aMtx2[0]+aMtx1[1]*aMtx2[1]+aMtx1[2]*aMtx2[2];
	tMtx[3] = aMtx1[0]*aMtx2[3]+aMtx1[1]*aMtx2[4]+aMtx1[2]*aMtx2[5];
	tMtx[6] = aMtx1[0]*aMtx2[6]+aMtx1[1]*aMtx2[7]+aMtx1[2]*aMtx2[8];
	tMtx[1] = aMtx1[3]*aMtx2[0]+aMtx1[4]*aMtx2[1]+aMtx1[5]*aMtx2[2];
	tMtx[4] = aMtx1[3]*aMtx2[3]+aMtx1[4]*aMtx2[4]+aMtx1[5]*aMtx2[5];
	tMtx[7] = aMtx1[3]*aMtx2[6]+aMtx1[4]*aMtx2[7]+aMtx1[5]*aMtx2[8];
	tMtx[2] = aMtx1[6]*aMtx2[0]+aMtx1[7]*aMtx2[1]+aMtx1[8]*aMtx2[2];
	tMtx[5] = aMtx1[6]*aMtx2[3]+aMtx1[7]*aMtx2[4]+aMtx1[8]*aMtx2[5];
	tMtx[8] = aMtx1[6]*aMtx2[6]+aMtx1[7]*aMtx2[7]+aMtx1[8]*aMtx2[8];
	
	double tVec6[3];
	tVec6[0] = tMtx[6]*cVec2[1] - tMtx[3]*cVec2[2];
	tVec6[1] = tMtx[0]*cVec2[2] - tMtx[6]*cVec2[0];
	tVec6[2] = tMtx[3]*cVec2[0] - tMtx[0]*cVec2[1];

	double tVec7[3];
	tVec7[0] = vVec1[0] - vVec2[0];
	tVec7[1] = vVec1[1] - vVec2[1];
	tVec7[2] = vVec1[2] - vVec2[2];

	double tVec8[3];
	tVec8[0] = aMtx1[0] * tVec7[0] + aMtx1[1] * tVec7[1] + aMtx1[2] * tVec7[2];
	tVec8[1] = aMtx1[3] * tVec7[0] + aMtx1[4] * tVec7[1] + aMtx1[5] * tVec7[2];
	tVec8[2] = aMtx1[6] * tVec7[0] + aMtx1[7] * tVec7[1] + aMtx1[8] * tVec7[2];

	double tVec9[3];
	CrossProduct(wVec1, cVec1, tVec9);

	VL[0] = tVec8[0] + tVec9[0];
	VL[1] = tVec8[1] +tVec9[1];
	VL[2] = tVec8[2] + tVec9[2];

	VL[0]-=wVec2[0]*tVec6[0];
	VL[0]-=wVec2[1]*tVec6[1];
	VL[0]-=wVec2[2]*tVec6[2];

	double tVec13[3];
	tVec13[0] = tMtx[7]*cVec2[1]-tMtx[4]*cVec2[2];
	tVec13[1] = tMtx[1]*cVec2[2]-tMtx[7]*cVec2[0];
	tVec13[2] = tMtx[4]*cVec2[0]-tMtx[1]*cVec2[1];

	VL[1]-=wVec2[0]*tVec13[0];
	VL[1]-=wVec2[1]*tVec13[1];
	VL[1]-=wVec2[2]*tVec13[2];

	double tVec14[3];
	tVec14[0] = tMtx[8]*cVec2[1]-tMtx[5]*cVec2[2];
	tVec14[1] = tMtx[2]*cVec2[2]-tMtx[8]*cVec2[0];
	tVec14[2] = tMtx[5]*cVec2[0]-tMtx[2]*cVec2[1];

	VL[2]-=wVec2[0]*tVec14[0];
	VL[2]-=wVec2[1]*tVec14[1];
	VL[2]-=wVec2[2]*tVec14[2];

	VL[3]=wVec1[0]-wVec2[0];
	VL[4]=wVec1[1]-wVec2[1];
	VL[5]=wVec1[2]-wVec2[2];


	//double tVec7[3];
	//tVec7[0]=-wVec1[2]*cVec1[1]+wVec1[1]*cVec1[2];
	//tVec7[1]=+wVec1[2]*cVec1[0]-wVec1[0]*cVec1[2];
	//tVec7[2]=-wVec1[1]*cVec1[0]+wVec1[0]*cVec1[1];

	//double tVec8[3];
	//tVec8[0]=-wVec1[2]*tVec7[1]+wVec1[1]*tVec7[2];
	//tVec8[1]=+wVec1[2]*tVec7[0]-wVec1[0]*tVec7[2];
	//tVec8[2]=-wVec1[1]*tVec7[0]+wVec1[0]*tVec7[1];

	//double tVec9[3];
	//tVec9[0]=-wVec2[2]*cVec2[1]+wVec2[1]*cVec2[2];
	//tVec9[1]=+wVec2[2]*cVec2[0]-wVec2[0]*cVec2[2];
	//tVec9[2]=-wVec2[1]*cVec2[0]+wVec2[0]*cVec2[1];

	//double tVec10[3];
	//tVec10[0]=-wVec2[2]*tVec9[1]+wVec2[1]*tVec9[2];
	//tVec10[1]=+wVec2[2]*tVec9[0]-wVec2[0]*tVec9[2];
	//tVec10[2]=-wVec2[1]*tVec9[0]+wVec2[0]*tVec9[1];

	//double tVec11[3];
	//tVec11[0]=tMtx[0]*tVec10[0]+tMtx[3]*tVec10[1]+tMtx[6]*tVec10[2];
	//tVec11[1]=tMtx[1]*tVec10[0]+tMtx[4]*tVec10[1]+tMtx[7]*tVec10[2];
	//tVec11[2]=tMtx[2]*tVec10[0]+tMtx[5]*tVec10[1]+tMtx[8]*tVec10[2];

	//// VEKMIN
	//double tVec12[3];
	//tVec12[0]=tVec11[0]-tVec8[0];
	//tVec12[1]=tVec11[1]-tVec8[1];
	//tVec12[2]=tVec11[2]-tVec8[2];
	//// VEKMIN END


	for (int i=3;i<6;i++)
		for (int nodeId1=0;nodeId1<6;nodeId1++)
		{
			A[i*6+nodeId1]=0.0;
			C[i*6+nodeId1]=0.0;
		}
//      a00 a01 a02 a03 a04 a05 
//      a06 a07 a08 a09 a10 a11
//      a12 a13 a14 a15 a16 a17
//      a18 a19 a20 a21 a22 a23
//      a24 a25 a26 a27 a28 a29
//      a30 a31 a32 a33 a34 a35 
//A
//       A0  A1  A2   0  X2 -X1 
//       A3  A4  A5 -X2   0  X0
//       A6  A7  A8  X1 -X0   0
//        0   0   0   1   0   0
//        0   0   0   0   1   0
//        0   0   0   0   0   1 
//C
// x=tVec6
// y=tVec13
// z=tVec14
//      -A0 -A1 -A2 -x1 -x2 -x3
//      -A3 -A4 -A5 -y1 -y2 -y3
//      -A6 -A7 -A8 -z1 -z2 -z3
//        0   0   0  -1   0   0
//        0   0   0   0  -1   0
//        0   0   0   0   0  -1 

	A[0]=aMtx1[0];
	A[1]=aMtx1[1];
	A[2]=aMtx1[2];
	A[3]=0.0;
	A[4]=cVec1[2];
	A[5]=-cVec1[1];
	A[6]=aMtx1[3];
	A[7]=aMtx1[4];
	A[8]=aMtx1[5];
	A[9]=-cVec1[2];
	A[10]=0.0;
	A[11]=cVec1[0];
	A[12]=aMtx1[6];
	A[13]=aMtx1[7];
	A[14]=aMtx1[8];
	A[15]=cVec1[1];
	A[16]=-cVec1[0];
	A[17]=0.0;
	A[21]=1.0;
	A[28]=1.0;
	A[35]=1.0;
	
	C[0]=-aMtx1[0];
    C[1]=-aMtx1[1];
    C[2]=-aMtx1[2];
	C[3]=-tVec6[0];
	C[4]=-tVec6[1];
	C[5]=-tVec6[2];
	C[6]=-aMtx1[3];
	C[7]=-aMtx1[4];
	C[8]=-aMtx1[5];
	C[9]=-tVec13[0];
	C[10]=-tVec13[1];
	C[11]=-tVec13[2];
	C[12]=-aMtx1[6];
	C[13]=-aMtx1[7];
	C[14]=-aMtx1[8];
	C[15]=-tVec14[0];
	C[16]=-tVec14[1];
	C[17]=-tVec14[2];
	C[21]=-1.0;
	C[28]=-1.0;
	C[35]=-1.0;
}

void StressStrainFortranSolver::CompareAngles
	(
		double* A1,
		double* A2,
		double* SL
	)
{
	double matrixAngleCosines[3];
	double linkAngleCosines[3];
	double angleDifferences[3];

	for (int i = 0; i < 3; i++)
	{
		matrixAngleCosines[i] = 0;
		for (int k = 0; k < 3; k++)
		{
			matrixAngleCosines[i] += A1[3 * i + k] * A2[3 * k + i];
		}
	}
	for (int i = 0; i < 3; i++)
	{
		linkAngleCosines[i] = cos(SL[i + 3]);
	}
	for (int i = 0; i < 3; i++)
	{
		angleDifferences[i] = fabs(matrixAngleCosines[i] - linkAngleCosines[i]);
	}
}

/**
* Получить векторный параметр
* @param coordinatesData - массив для записи векторного параметра
* @param scaleFactor - масштабный коэффициент
*/
void StressStrainFortranSolver::GetVectorParameter
	(
		float* coordinatesData,
		double scaleFactor
	)	const
{
	for (int nodeIndex = 0; nodeIndex < _nElements; ++nodeIndex)
	{
		const double* initialNodeCoordinates = _nodes + 3 * nodeIndex;
		const double* currentNodeCoordinates = _dataInternal + 6 * nodeIndex;
		float* nodeCoordinatesData = coordinatesData + 3 * nodeIndex;

		for (int directionIndex = 0; directionIndex < 3; ++directionIndex)
		{
			const double initialNodeCoordinate = initialNodeCoordinates[directionIndex];
			const double currentNodeCoordinate = currentNodeCoordinates[directionIndex];
			float& nodeCoordinateData = nodeCoordinatesData[directionIndex];

			if (fabs(currentNodeCoordinate) > 1e20)
			{
				std::cout << "COORDINATE BIG NUMBER!!!" << std::endl;
			}
			if (scaleFactor == 1.0)
			{
				nodeCoordinateData = static_cast<float>(currentNodeCoordinate);
			}
			else
			{
				nodeCoordinateData =
					static_cast<float>
						(
							initialNodeCoordinate +
							(currentNodeCoordinate - initialNodeCoordinate) * scaleFactor
						);
			}
		}
	}
}