#include "StressStrainCppSolver.h"

#include "CsrSymmetricMatrix.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"

#include <omp.h>
#include "../../AdditionalModules/fmath/Matrix3x4.h"


using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;


//#define NO_LINKSH
//#define NO_INTOMSUB


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


namespace Stress
{

	void PrintVector(Vec3 vec, string comment)
	{
		std::cout << comment << "{";
		for (int i = 0; i < 3; i++)
			std::cout << vec[i] << ' ';
		std::cout << "}\n";
	}

	void PrintMatrix(Mat3 mat, string comment)
	{
		std::cout << comment << std::endl;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				std::cout << mat[i * 3 + j] << ' ';
			std::cout << std::endl;
		}
		//std::cout << "}\n";
	}

	void PrintMatrix(MathHelpers::Mat3x4 mat, string comment)
	{
		std::cout << comment << std::endl;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				std::cout << mat.E(i,j) << ' ';
			std::cout << std::endl;
		}
		//std::cout << "}\n";
	}


#define MAX(x, y) ((x) > (y) ? x : y)


double DMOD
	(
		const double d1,
		const double d2
	)
{
	return d1 - d2 * ((int)(d1 / d2));
}

StressStrainCppSolver::StressStrainCppSolver
	(
		double* params,		// параметры
		int* links,			// индексы элементов связей (пары)
		int nLinks,			// число связей
		double *nodes,		// координаты узлов
		int nNodes,			// число узлов
		double gridStep,	// шаг сетки
		double timeStep,	// шаг по времени
		int numThreads		// число потоков
	)
	:
		StressStrainSolver(nNodes),
		_isFirstSolution(true),
		_isFirstIteration(true),
		_nIteration(0),
		_nElements(nNodes),
		_poissonRatio(0.),
		_stageRK(0),                // этап вычисления итерации РК4 (0-3)
		_time(0),					// время
		_timeStep(timeStep),		// шаг интегрирования
		_timeStep2(timeStep*0.5),   // половина шага интегрирования
		_timeStep4(timeStep*0.25),  // четверть шага интегрирования
		_gridStep(gridStep),		// шаг сетки
		_elasticModulus(params[0]),	// модуль упругости
		_dampingFactor(params[1]),	// коэффициент демпфирования
		_density(params[2]),		// плотность материала
		_stiffScale(params[3])

{
	std::cout << "------------- CPP SOLVER IS CREATED ---------------------------" << std::endl
		<< "VECSTRIDE =" << vecStride << std::endl << std::endl;
	
	_nVariables = nNodes * vecStride2;

	// масштабный коэффициент
	_gridStep2 = _gridStep * _gridStep;
	_gridStep3 = _gridStep * _gridStep2;

	_cellMass = _density * _gridStep3; // масса ячейки
	_elasticModulusScaled = _elasticModulus / _stiffScale; // отмасштабированный модуль упругости

	_numThreads = (numThreads > 0) ? numThreads : omp_get_max_threads();
	
	_shearModulusScaled = _elasticModulusScaled / 2;
	_dampingFactorLinear = 0.01 * 2 * _dampingFactor * sqrt(_elasticModulusScaled * _cellMass * _gridStep);
	_dampingFactorAngular = 10 * _dampingFactor * sqrt(2 * _elasticModulusScaled * _cellMass * _gridStep / 3) * _gridStep2;
	
	std::cout
		<< "################### STIFF SCALE: " << _stiffScale
		<< " ################### NUMBER OF THREADS: " << _numThreads << " param=" << numThreads  << std::endl;

// Выделение памяти

	size_t varSize = _nVariables * sizeof(double);


	_initX = (double*)_aligned_malloc(varSize,alignment);
	_initDX= (double*)_aligned_malloc(varSize,alignment);
	_hDDX1 = (double*)_aligned_malloc(varSize,alignment);
	_hDDX2 = (double*)_aligned_malloc(varSize,alignment);
	_hDDX3 = (double*)_aligned_malloc(varSize,alignment);

	_varX = _dataInternal;
	_varDX = _dataInternal+_nVariables;
	_varDDX = _dataInternal+_nVariables*2;

	//R = new double[_nVariables];
	//RZ = new double[_nVariables];
	//R1Z = new double[_nVariables];
	_fue = new Fue(nNodes);
	_copym = new Copym(nNodes);

	// Связанные элементы: порядок x-,y-,z-,x+,y+,z+ 
	// Нумерация с 1, если 0, то связь отсутствует
	//    1
	//    |
	// 0 -x- 3
	//    |
	//    4
	// [2] = z-, [5] = z+
	_linkedElements = new int[nNodes * 6];

	// копируем координаты узлов сетки
	_elements = new double[nNodes*3];
	memcpy(_elements, nodes, sizeof(double) * nNodes * 3);
	memset(_linkedElements, 0, nNodes * 6 * sizeof(int));
	for (int i = 0; i < nLinks * 2; i += 2)
	{
		// номера элементов в связях
		int elementId1 = links[i];
		int elementId2 = links[i + 1];
		int* linkedElementsOffset1 = _linkedElements + elementId1 * 6;
		int* linkedElementsOffset2 = _linkedElements + elementId2 * 6;

		// выбираем одну из координат
		double* node1 = &_elements[3 * links[i]];
		double* node2 = &_elements[3 * links[i + 1]];

		const double relativeTolerance = 1e-4;
		for (int k = 0; k < 3; k++)
		{
			double dx = node1[k] - node2[k];
			if (dx < 0 && fabs(gridStep + dx) < gridStep * relativeTolerance)
			{
				linkedElementsOffset2[k] = elementId1 + 1;
				linkedElementsOffset1[k + 3] = elementId2 + 1;
			}
			else if (dx > 0 && fabs(gridStep - dx) < gridStep * relativeTolerance)
			{
				linkedElementsOffset2[k + 3] = elementId1 + 1;
				linkedElementsOffset1[k] = elementId2 + 1;
			}
		}
	}
	memset(_dataInternal, 0, _nVariables * 3 * sizeof(double));
	for (int i = 0; i < nNodes; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_dataInternal[i * 2 * vecStride + j] = nodes[i * 3 + j];
		}
	}
	intomsub();
	_testTimer.Allocate(10);
}

// virtual
StressStrainCppSolver::~StressStrainCppSolver()
{
	//delete [] GR1;
	//_aligned_free(_rotationMatrices);
	_aligned_free(_initX);
	_aligned_free(_initDX);
	_aligned_free(_hDDX1);
	_aligned_free(_hDDX2);
	_aligned_free(_hDDX3);

	//delete [] R;
	//delete [] RZ;
	//delete [] R1Z;
	delete [] _linkedElements;
	delete [] _elements;
	delete _fue;
	delete _copym;
}

/** Добавление связей
* @param links - массив 6*_nElements с признаками связей (1 - есть связь, 0 - нет связи)
*/
// virtual
void StressStrainCppSolver::AddLinks
	(
		const int* links
	)
{
	for (int i = 0; i < _nElements * 6; i++)
	{
		_linkedElements[i] = links[i];
	}
}

void StressStrainCppSolver::AddPartialBoundary
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
void StressStrainCppSolver::AddBoundary
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
			numberOfBoundaryNodes,
			vecStride
		);

	_boundaryParamsSet.push_back(bp);
}

// virtual
void StressStrainCppSolver::ChangeBoundaryParam
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
double StressStrainCppSolver::GetBoundaryParam
	(
		const int bcNumber,
		const int bcParamNumber
	)
{
	return _boundaryParamsSet[bcNumber].GetParam(bcParamNumber);
}

// virtual
void StressStrainCppSolver::UpdateBuffer
	(
		double scale
	)
{
	int id = 0;

	GetScalarParameter(_data);
	for (int i = 0; i < _nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (fabs(_dataInternal[i * vecStride2 + j]) > 1e20)
			{
				std::cout << "COORDINATE BIG NUMBER!!!" << std::endl;
			}
			_data[i * 3 + j + _nElements] =
				(float)
					(
						_elements[i * 3 + j] +
						(_dataInternal[i * vecStride2 + j] - _elements[i * 3 + j]) * scale
					);
		}
	}
}

void StressStrainCppSolver::intomsub()
{
#ifdef NOINTOMSUB
	return;
#endif

//#pragma omp parallel for private(i) num_threads(NumThreads)
	//std::cout << "                _nElements = " << _nElements << std::endl;
	for (int i = 0; i < _nElements; i++)
	{
		urejl4s
			(
				&_dataRotationMtx[i * matStride],
				&_dataInternal[i * vecStride2 + vecStride],
				&_dataInternal[i * vecStride2 + vecStride + _nElements * vecStride2],
				_stageRK,
				i
			);
	}
}

void StressStrainCppSolver::urejl4s
	(
		double* a,
		double *ug,
		double *om,
		int mets,
		const int id
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	int K = id * vecStride;
	int KA = id * matStride;
	
	if (mets == 0) 
	{
		_fue->Update(K);
		_fue->FLAGRY[id] = 0;
		_copym->Copy(a, _fue->A1 + KA, id, _time);
		
		//std::cout << "#############################################" << std::endl;
		//std::cout << "intomsub urejl4s METS == 0" << std::endl;

		return;
	}
	
	double UY = _fue->UE[K];
	UY = fabs(DMOD(UY, 2 * M_PI));           // DMOD - ???????

	if ((mets == 1) && (UY > 1.28))
	{

		_fue->Update(K);
		_fue->FLAGRY[id] = 1.0;
		_fue->IFLAGRY = 1;
		
		_copym->Copy(a, _fue->A1 + KA, id, _time);
	} 
	_fue->Update(K, mets);
	_fue->UpdateR(K, om, _timeStep);
	_fue->UpdateR2(K, mets);

	_fue->UpdateMtx(K, a);
	MatrixMul(_fue->A1 + KA, a);
}

void StressStrainCppSolver::MatrixMul
	(
		double *a1,
		double *a2
	)
{
	double tmp[matStride];
	memcpy(&tmp[0], a2, sizeof(double)*matStride);
	for (int row = 0; row < 3; row++)
	{
		for (int col = 0; col < 3; col++)  
		{
			a2[row * vecStride + col] = 0.;
			for(int k = 0; k < 3; k ++)
			{
				a2[row * vecStride + col] += a1[row * vecStride + k] * tmp[vecStride * k + col];
			}
		}
	}

	//Mat3 mtx1(a1);
	//Mat3 mtx2(a2);
	//Mat3 mtx3 = mtx1 * mtx2;
	//mtx3.Export(a2);
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

void StressStrainCppSolver::linksh4AVX
(
size_t side,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *SL,		// выход деформаций
double *VL,		// выход изм. скоростей
double& rx,		// выход
double& ry,		// выход 
double& rz,		// выход 
int nodeId1,	// номер узла 1
int nodeId2,	// номер узла 2
int nNodes		// количество узлов
)
{
	__declspec(align(32)) double nodeVectors[] =
	{
		-1, 0, 0, 0,
		 0,-1, 0, 0,
		 0, 0,-1, 0,
		 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0, 1, 0
	};

	double* pNodeVectors = &nodeVectors[0] + side;
	for (size_t component = 0; component < 3; component++)
		pNodeVectors[component] *= (_gridStep * 0.5);
	rx = pNodeVectors[0];
	ry = pNodeVectors[1];
	rz = pNodeVectors[2];

	size_t strideDeriv = nNodes * vecStride2;
	size_t nodeOffset1 = nodeId1 * vecStride2;
	size_t nodeOffset2 = nodeId2 * vecStride2;

	// Start AVX code
	double* pmatA01 = _dataRotationMtx + nodeId1 * matStride;
	double* pmatA02 = _dataRotationMtx + nodeId2 * matStride;

	// vecStride must be 4
	__m256d matA01row1 = _mm256_load_pd(pmatA01);
	__m256d matA01row2 = _mm256_load_pd(pmatA01 + vecStride);
	__m256d matA01row3 = _mm256_load_pd(pmatA01 + vecStride2);

	__m256d matA02el1;
	__m256d matA02el2;
	__m256d matA02el3;

	__m256d matA21row1;
	__m256d matA21row2;
	__m256d matA21row3;

	// matA02 column 1/3
	matA02el1 = _mm256_set1_pd(pmatA02[0]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2]);

	// matA01.TMul(matA02)
	matA02el1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA01row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA01row3, matA02el3);

	matA21row1 = _mm256_add_pd(matA02el1, _mm256_add_pd(matA02el2, matA02el3));

	// matA02 column 2/3
	matA02el1 = _mm256_set1_pd(pmatA02[1]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 1]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 1]);

	// matA01.TMul(matA02)
	matA02el1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA01row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA01row3, matA02el3);

	matA21row2 = _mm256_add_pd(matA02el1, _mm256_add_pd(matA02el2, matA02el3));

	// matA02 column 3/3
	matA02el1 = _mm256_set1_pd(pmatA02[2]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 2]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 2]);

	// matA01.TMul(matA02)
	matA02el1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA01row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA01row3, matA02el3);

	matA21row3 = _mm256_add_pd(matA02el1, _mm256_add_pd(matA02el2, matA02el3));
	// Матрица A_{21} сформирована

#ifdef DEBUG_OUT
	// выгрузка для отладки
	__declspec(align(32)) double matrix[12];
	_mm256_store_pd(&matrix[0], matA21row1);
	_mm256_store_pd(&matrix[4], matA21row2);
	_mm256_store_pd(&matrix[8], matA21row3);

	Mat3x4 mtx(matrix);
	PrintMatrix(mtx, "AVX A21");
#endif

	__m256d ivecC1 = _mm256_load_pd(pNodeVectors);
	__m256d vecDP = _mm256_sub_pd(
		_mm256_load_pd(_dataInternal + nodeOffset1),
		_mm256_load_pd(_dataInternal + nodeOffset2)); // P1-P2

	matA02el1 = _mm256_set1_pd(pNodeVectors[0]);
	matA02el2 = _mm256_set1_pd(pNodeVectors[1]);
	matA02el3 = _mm256_set1_pd(pNodeVectors[2]);

	matA02el1 = _mm256_mul_pd(matA21row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA21row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA21row3, matA02el3);

	__m256d mul1 = _mm256_add_pd(_mm256_add_pd(matA02el1, matA02el2), matA02el3);

	matA02el1 = _mm256_set1_pd(vecDP.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(vecDP.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(vecDP.m256d_f64[2]);

	matA02el1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA01row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA01row3, matA02el3);


	__m256d mul2 = _mm256_add_pd(_mm256_add_pd(matA02el1, matA02el2), matA02el3);
	__m256d res = _mm256_add_pd(_mm256_add_pd(ivecC1, mul1), mul2);

	_mm256_storeu_pd(SL, res); // получено SL, линейные компоненты
	
	// Расчет VL
	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecC1.Cross(vecW1) - matA12*(vecC2.Cross(vecW1));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	__declspec(align(32)) double cp1[4] = {0};	// векторное произведение

	CrossProduct(pNodeVectors, _dataInternal + nodeOffset1 + strideDeriv + 3, &cp1[0]);	// [cVec x w]
	__m256d cp1r = _mm256_load_pd(&cp1[0]);
	__m256d vecDV = _mm256_sub_pd(
		_mm256_load_pd(_dataInternal + nodeOffset1 + strideDeriv),
		_mm256_load_pd(_dataInternal + nodeOffset2 + strideDeriv)); // V1-V2

	matA02el1 = _mm256_set1_pd(cp1[0]);
	matA02el2 = _mm256_set1_pd(cp1[1]);
	matA02el3 = _mm256_set1_pd(cp1[2]);

	matA02el1 = _mm256_mul_pd(matA21row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA21row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA21row3, matA02el3);

	mul1 = _mm256_add_pd(_mm256_add_pd(matA02el1, matA02el2), matA02el3);

	matA02el1 = _mm256_set1_pd(vecDV.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(vecDV.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(vecDV.m256d_f64[2]);

	matA02el1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA02el2 = _mm256_mul_pd(matA01row2, matA02el2);
	matA02el3 = _mm256_mul_pd(matA01row3, matA02el3);

	mul2 = _mm256_add_pd(_mm256_add_pd(matA02el1, matA02el2), matA02el3);

	res = _mm256_add_pd(_mm256_add_pd(cp1r, mul1), mul2);

	_mm256_storeu_pd(VL, res); // получено VL, линейные компоненты

	double *pVec1 = _dataInternal + nodeOffset1;
	double *pVec2 = _dataInternal + nodeOffset2;
	Vec3Ref vecR1 = MakeVec3(pVec1 + vecStride);
	Vec3Ref vecR2 = MakeVec3(pVec2 + vecStride);
	Vec3Ref vecW1 = MakeVec3(pVec1 + strideDeriv + 3);
	Vec3Ref vecW2 = MakeVec3(pVec2 + strideDeriv + 3);
	//TODO: для угловых степеней свободы - просто берутся разница углов и угловых скоростей
	(vecR1 - vecR2).Export(SL + vecStride);
	(vecW1 - vecW2).Export(VL + vecStride);
}

void StressStrainCppSolver::linksh4
(
size_t side,	// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *SL,		// выход деформаций
double *VL,		// выход изм. скоростей
double& rx,		// выход
double& ry,		// выход 
double& rz,		// выход 
int nodeId1,	// номер узла 1
int nodeId2,	// номер узла 2
int nNodes		// количество узлов
)
{
	__declspec(align(32)) double nodeVectors[] =
	{	
	   -1, 0, 0, 0,
		0,-1, 0, 0,
		0, 0,-1, 0,
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0
	};

	double* pNodeVectors = &nodeVectors[0] + side;

	//unsigned int strideVec = 3;
	//unsigned int strideMat = strideVec * strideVec;
	size_t strideDeriv = nNodes * vecStride2;
	size_t nodeOffset1 = nodeId1 * vecStride2;
	size_t nodeOffset2 = nodeId2 * vecStride2;

	MathHelpers::Mat3x4 matA01(_dataRotationMtx + nodeId1 * matStride);
	MathHelpers::Mat3x4 matA02(_dataRotationMtx + nodeId2 * matStride);
	//Mat3 matA12 = matA01.Tmul(matA02);
	MathHelpers::Mat3x4 matA21 = matA02.Tmul(matA01);

	Vec3 vecC1 = MakeVec3(pNodeVectors)*_gridStep * 0.5;
	Vec3 vecC2 = -vecC1;// MakeVec3(cVec2);

	rx = vecC1[0];
	ry = vecC1[1];
	rz = vecC1[2];

	double *pVec1 = _dataInternal + nodeOffset1;
	double *pVec2 = _dataInternal + nodeOffset2;
	Vec3Ref vecP1 = MakeVec3(pVec1);
	Vec3Ref vecP2 = MakeVec3(pVec2);
	Vec3Ref vecR1 = MakeVec3(pVec1 + vecStride);
	Vec3Ref vecR2 = MakeVec3(pVec2 + vecStride);
	Vec3Ref vecV1 = MakeVec3(pVec1 + strideDeriv);
	Vec3Ref vecV2 = MakeVec3(pVec2 + strideDeriv);
	Vec3Ref vecW1 = MakeVec3(pVec1 + strideDeriv + 3);
	Vec3Ref vecW2 = MakeVec3(pVec2 + strideDeriv + 3);
#ifdef DEBUG_OUT
	PrintMatrix(matA01, "A01=");
	PrintMatrix(matA02, "A02=");
	PrintMatrix(matA21, "A21=");
#endif
	// переводим вектор линии точек связи С2-С1 в СК1
	Vec3 vecT0 = vecC1 - matA21.Tmul(vecC2) - matA01.Tmul(vecP2 - vecP1);

	Vec3 v1 = matA21.Tmul(vecC2);
	Vec3 v2 = matA01.Tmul(vecP2 - vecP1);

	vecT0.Export(SL);
	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecC1.Cross(vecW1) - matA21.Tmul(vecC2.Cross(vecW1));
	VecT1.Export(VL);

	Vec3 v3 = vecC2.Cross(vecW1);
	Vec3 v4 = matA21.Tmul(v3);

#ifdef DEBUG_OUT
	PrintVector(v1, "v1=");
	PrintVector(v2, "v2=");
	PrintVector(v3, "v3=");
	PrintVector(v4, "v4=");
#endif



	// для угловых степеней свободы - просто берутся разница углов и угловых скоростей
	(vecR1 - vecR2).Export(SL + vecStride);
	(vecW1 - vecW2).Export(VL + vecStride);


	//Vec3 vecT2 = matA12.Rcross(0, vecC2);
	//Vec3 vecT3 = matA12.Rcross(1, vecC2);
	//Vec3 vecT4 = matA12.Rcross(2, vecC2);

	//vecT2.Export(C + 3);
	//vecT3.Export(C + 9);
	//vecT4.Export(C + 15);

	// Заполнение матриц A и C
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
	// x=A12.Row[0] x VecC2
	// y=A12.Row[1] x VecC2
	// z=A12.Row[2] x VecC2
	//      -A0 -A1 -A2 -x1 -x2 -x3
	//      -A3 -A4 -A5 -y1 -y2 -y3
	//      -A6 -A7 -A8 -z1 -z2 -z3
	//        0   0   0  -1   0   0
	//        0   0   0   0  -1   0
	//        0   0   0   0   0  -1 

	//for (int i = 3; i<6; i++)
	//	for (int nodeId1 = 0; nodeId1<6; nodeId1++)
	//	{
	//		A[i * 6 + nodeId1] = 0.0;
	//		C[i * 6 + nodeId1] = 0.0;
	//	}

	//for (int i = 0; i < 3; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		A[i * 6 + j] = matA01[i * 3 + j];
	//		C[i * 6 + j] = -matA01[i * 3 + j];
	//	}
	//}

	//A[3] = 0.0;        A[4] = vecC1[2];   A[5] = -vecC1[1];
	//A[9] = -vecC1[2];  A[10] = 0.0;       A[11] = vecC1[0];
	//A[15] = vecC1[1];  A[16] = -vecC1[0]; A[17] = 0.0;

	//A[21] = 1.0;
	//A[28] = 1.0;
	//A[35] = 1.0;
}


void StressStrainCppSolver::linksh3
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
	unsigned int strideVec = 3;
	unsigned int strideMat = strideVec * strideVec;
	unsigned int strideDeriv = nNodes * 2 * strideVec;
	unsigned int nodeOffset1 = nodeId1 * 2 * strideVec;
	unsigned int nodeOffset2 = nodeId2 * 2 * strideVec;
	
	Mat3 matA01(_dataRotationMtx + nodeId1 * strideMat);
	Mat3 matA02(_dataRotationMtx + nodeId2 * strideMat);
	Vec3Ref vecC1 = MakeVec3(cVec1);
	Vec3Ref vecC2 = MakeVec3(cVec2);
	Mat3 matA12 = matA01.Tmul(matA02);

	double *pVec1 = _dataInternal + nodeOffset1;
	double *pVec2 = _dataInternal + nodeOffset2;
	Vec3Ref vecP1 = MakeVec3(pVec1);
	Vec3Ref vecP2 = MakeVec3(pVec2);
	Vec3Ref vecR1 = MakeVec3(pVec1 + strideVec);
	Vec3Ref vecR2 = MakeVec3(pVec2 + strideVec);
	Vec3Ref vecV1 = MakeVec3(pVec1 + strideDeriv);
	Vec3Ref vecV2 = MakeVec3(pVec2 + strideDeriv);
	Vec3Ref vecW1 = MakeVec3(pVec1 + strideDeriv + 3);
	Vec3Ref vecW2 = MakeVec3(pVec2 + strideDeriv + 3);

	// переводим вектор линии точек связи С2-С1 в СК1
	Vec3 vecT0 = vecC1 - matA12*vecC2 - matA01.Tmul(vecP2-vecP1);
	vecT0.Export(SL);

	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	Vec3 VecT1 = matA01.Tmul(vecV1-vecV2) + vecC1.Cross(vecW1) - matA12*(vecC2.Cross(vecW1));
	VecT1.Export(VL);

	// для угловых степеней свободы - просто берутся разница углов и угловых скоростей
	(vecR1-vecR2).Export(SL+3);
	(vecW1-vecW2).Export(VL+3);


	Vec3 vecT2 = matA12.Rcross(0, vecC2);
	Vec3 vecT3 = matA12.Rcross(1, vecC2);
	Vec3 vecT4 = matA12.Rcross(2, vecC2);

	vecT2.Export(C+3);
	vecT3.Export(C+9);
	vecT4.Export(C+15);

	// Заполнение матриц A и C
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
	// x=A12.Row[0] x VecC2
	// y=A12.Row[1] x VecC2
	// z=A12.Row[2] x VecC2
	//      -A0 -A1 -A2 -x1 -x2 -x3
	//      -A3 -A4 -A5 -y1 -y2 -y3
	//      -A6 -A7 -A8 -z1 -z2 -z3
	//        0   0   0  -1   0   0
	//        0   0   0   0  -1   0
	//        0   0   0   0   0  -1 

	for (int i=3;i<6;i++)
		for (int nodeId1=0;nodeId1<6;nodeId1++)
		{
			A[i*6+nodeId1]=0.0;
			C[i*6+nodeId1]=0.0;
		}

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			A[i*6+j] = matA01[i*3+j];
			C[i*6+j] = -matA01[i*3+j];
		}
	}
	
	A[3]=0.0;        A[4]=vecC1[2];   A[5]=-vecC1[1]; 
	A[9]=-vecC1[2];  A[10]=0.0;       A[11]=vecC1[0];
	A[15]=vecC1[1];  A[16]=-vecC1[0]; A[17]=0.0;

	A[21]=1.0;
	A[28]=1.0;
	A[35]=1.0;
}

void StressStrainCppSolver::linksh2
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

	double *aMtx1 = _dataRotationMtx + nodeId1 * 9;
	double *aMtx2 = _dataRotationMtx + nodeId2 * 9;

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

	VL[0]=vVec1[0]*aMtx1[0]+vVec1[1]*aMtx1[1]+vVec1[2]*aMtx1[2];
	VL[0]+=wVec1[1]*cVec1[2]-wVec1[2]*cVec1[1];
	VL[0]-=vVec2[0]*aMtx1[0];
	VL[0]-=vVec2[1]*aMtx1[1];
	VL[0]-=vVec2[2]*aMtx1[2];
	VL[0]-=wVec2[0]*tVec6[0];
	VL[0]-=wVec2[1]*tVec6[1];
	VL[0]-=wVec2[2]*tVec6[2];

	double tVec13[3];
	tVec13[0] = tMtx[7]*cVec2[1]-tMtx[4]*cVec2[2];
	tVec13[1] = tMtx[1]*cVec2[2]-tMtx[7]*cVec2[0];
	tVec13[2] = tMtx[4]*cVec2[0]-tMtx[1]*cVec2[1];

	VL[1]=vVec1[0]*aMtx1[3]+vVec1[1]*aMtx1[4]+vVec1[2]*aMtx1[5];
	VL[1]+=wVec1[2]*cVec1[0]-wVec1[0]*cVec1[2];
	VL[1]-=vVec2[0]*aMtx1[3];
	VL[1]-=vVec2[1]*aMtx1[4];
	VL[1]-=vVec2[2]*aMtx1[5];
	VL[1]-=wVec2[0]*tVec13[0];
	VL[1]-=wVec2[1]*tVec13[1];
	VL[1]-=wVec2[2]*tVec13[2];

	double tVec14[3];
	tVec14[0] = tMtx[8]*cVec2[1]-tMtx[5]*cVec2[2];
	tVec14[1] = tMtx[2]*cVec2[2]-tMtx[8]*cVec2[0];
	tVec14[2] = tMtx[5]*cVec2[0]-tMtx[2]*cVec2[1];

	VL[2]=vVec1[0]*aMtx1[6]+vVec1[1]*aMtx1[7]+vVec1[2]*aMtx1[8];
	VL[2]+=wVec1[0]*cVec1[1]-wVec1[1]*cVec1[0];
	VL[2]-=vVec2[0]*aMtx1[6];
	VL[2]-=vVec2[1]*aMtx1[7];
	VL[2]-=vVec2[2]*aMtx1[8];
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



void StressStrainCppSolver::FindStressStrainMatrix()
{
	double lambda = _poissonRatio / (1. - 2*_poissonRatio)/(1. + _poissonRatio);
	double mu = 0.5 / (1. + _poissonRatio);

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			_lameMatrix[i][j] = 0;
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_lameMatrix[i][j] = lambda;
			if (i == j)
			{
				_lameMatrix[i][j] += mu * 2.;
			}
		}
	}
	for (int i = 3; i < 6; i++)
	{
		_lameMatrix[i][i] = mu;
	}
}
}