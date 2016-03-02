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


StressStrainCppSolver::StressStrainCppSolver
	(
		double* params,		// ���������
		int* links,			// ������� ��������� ������ (����)
		int nLinks,			// ����� ������
		double *nodes,		// ���������� �����
		int nElements,		// ����� �����
		double gridStep,	// ��� �����
		double timeStep,	// ��� �� �������
		int numThreads,		// ����� �������
		int stride
	)
	:
		StressStrainSolver(nElements, stride),
		_isFirstSolution(true),
		_isFirstIteration(true),
		_nIteration(0),
		_poissonRatio(0.),
		_stageRK(0),                // ���� ���������� �������� ��4 (0-3)
		_time(0),					// �����
		_timeStep(timeStep),		// ��� ��������������
		_timeStep2(timeStep*0.5),   // �������� ���� ��������������
		_timeStep4(timeStep*0.25),  // �������� ���� ��������������
		_gridStep(gridStep),		// ��� �����
		_elasticModulus(params[0]),	// ������ ���������
		_dampingFactor(params[1]),	// ����������� �������������
		_density(params[2]),		// ��������� ���������
		_stiffScale(params[3])

{
	
	_nVariables = nElements * vecStride2;

	// ���������� �����������
	_gridStep2 = _gridStep * _gridStep;
	_gridStep3 = _gridStep * _gridStep2;

	_cellMass = _density * _gridStep3; // ����� ������
	_elasticModulusScaled = _elasticModulus / _stiffScale; // ������������������ ������ ���������

	_numThreads = (numThreads > 0) ? numThreads : omp_get_max_threads();
	
	_shearModulusScaled = _elasticModulusScaled / 2;
	_dampingFactorLinear = 0.01 * 2 * _dampingFactor * sqrt(_elasticModulusScaled * _cellMass * _gridStep);
	_dampingFactorAngular = 10 * _dampingFactor * sqrt(2 * _elasticModulusScaled * _cellMass * _gridStep / 3) * _gridStep2;
	
	const int outWidth = 15;
	std::cout << "------------------------------" << std::endl
			  << "    CPP SOLVER IS CREATED" << std::endl
			  << "------------------------------" << std::endl
		<< std::setw(outWidth) << "THREADS: " << std::setw(outWidth) << _numThreads << std::endl
		<< std::setw(outWidth) << "VECSTRIDE: " << std::setw(outWidth) << vecStride << std::endl
		<< std::setw(outWidth) << "STIFF SCALE: " << std::setw(outWidth) << _stiffScale << std::endl
		<< std::setw(outWidth) << "ELASTIC: "   << std::setw(outWidth) << _elasticModulus << std::endl
		<< std::setw(outWidth) << "DAMPING: "   << std::setw(outWidth) << _dampingFactor << std::endl
		<< std::setw(outWidth) << "DENSITY: "   << std::setw(outWidth) << _density << std::endl
		<< std::setw(outWidth) << "ELEMENTS: "  << std::setw(outWidth) << nElements << std::endl
		<< std::setw(outWidth) << "GRID STEP: " << std::setw(outWidth) << _gridStep << std::endl;

	
	// ��������� ������

	size_t varSize = _nVariables * sizeof(double);


	_initX = (double*)_aligned_malloc(varSize,alignment);
	_initDX= (double*)_aligned_malloc(varSize,alignment);
	_hDDX1 = (double*)_aligned_malloc(varSize,alignment);
	_hDDX2 = (double*)_aligned_malloc(varSize,alignment);
	_hDDX3 = (double*)_aligned_malloc(varSize,alignment);



	_varX = _dataInternal;
	_varDX = _dataInternal+_nVariables;
	_varDDX = _dataInternal+_nVariables*2;

	//_varDDX = new double[_nVariables];
	//RZ = new double[_nVariables];
	//R1Z = new double[_nVariables];
	_rotationSolver = new RotationSolver(nElements, stride);
	_copym = new Copym(nElements, stride);

	// ��������� ��������: ������� x-,y-,z-,x+,y+,z+ 
	// ��������� � 1, ���� 0, �� ����� �����������
	//    1
	//    |
	// 0 -x- 3
	//    |
	//    4
	// [2] = z-, [5] = z+
	_linkedElements = new int[nElements * 6];
	_radiusVectors = (double*)_aligned_malloc(vecStride * 6 * sizeof(double), alignment);
	memset(_radiusVectors, 0, vecStride * 6 * sizeof(double));

	for (int i = 0; i < 3; i++)
	{
		_radiusVectors[i*vecStride*2 + i] = 0.5 * _gridStep;
		_radiusVectors[i*vecStride + i] = -0.5 * _gridStep;
	}

	// �������� ���������� ����� �����
	_elements = new double[nElements*3];
	memcpy(_elements, nodes, sizeof(double) * nElements * 3);
	memset(_linkedElements, 0, nElements * 6 * sizeof(int));
	for (int i = 0; i < nLinks * 2; i += 2)
	{
		// ������ ��������� � ������
		int elementId1 = links[i];
		int elementId2 = links[i + 1];
		int* linkedElementsOffset1 = _linkedElements + elementId1 * 6;
		int* linkedElementsOffset2 = _linkedElements + elementId2 * 6;

		// �������� ���� �� ���������
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
	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_dataInternal[i * 2 * vecStride + j] = nodes[i * 3 + j];
		}
	}
	CalculateRotations();
	_testTimer.Allocate(10);
}

// virtual
StressStrainCppSolver::~StressStrainCppSolver()
{
	_aligned_free(_initX);
	_aligned_free(_initDX);
	_aligned_free(_hDDX1);
	_aligned_free(_hDDX2);
	_aligned_free(_hDDX3);
	_aligned_free(_radiusVectors);

	delete [] _linkedElements;
	delete [] _elements;
	delete _rotationSolver;
	delete _copym;
}

/** ���������� ������
* @param links - ������ 6*_nElements � ���������� ������ (1 - ���� �����, 0 - ��� �����)
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
		// TODO: �������� exception
		;	// ������!!!
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

void StressStrainCppSolver::CalculateRotations()
{
#ifdef NOINTOMSUB
	return;
#endif

//#pragma omp parallel for private(i) num_threads(NumThreads)
	//std::cout << "                _nElements = " << _nElements << std::endl;
	for (size_t elementId = 0; elementId < _nElements; elementId++)
	{
		SolveElementRotation
		(
			GetRotationMatrix(elementId),
			GetElementVelocityAngular(elementId),
			elementId
		);
	}
}

void StressStrainCppSolver::SolveElementRotation
	(
		double* elementMtx,
		double *elementW,
		int elementId
	)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	// �������� �� ����������� �������
	double* rframeMtx = _rotationSolver->GetRframeMtx(elementId);

	// ������������� �� ������� ����
	if (_stageRK == 0) 
	{
		_rotationSolver->MakeZeroVectors(elementId);
		_copym->Copy(elementMtx, rframeMtx, elementId, _time);		
		return;
	}
	
	// �������� ������������� �� ������ ����
	if ((_stageRK == 1) && _rotationSolver->IsSingularityAngle(elementId))
	{
		_rotationSolver->MakeZeroVectors(elementId);
		_copym->Copy(elementMtx, rframeMtx, elementId, _time);
		// TODO: memcpy(rframeMtx, elementMtx, sizeof(double) * matStride);
	} 
	
	// ���������� �����
	_rotationSolver->Update(elementId, elementW, elementMtx, _timeStep, _stageRK);
	
	// ��������� � ����� �� ��� �������������
	MatrixMul(rframeMtx, elementMtx);
}

void StressStrainCppSolver::MatrixMul
	(
		double *a1,
		double *a2
	)
{
	if (vecStride == 3)
	{
		Mat3 mtx1(a1);
		Mat3 mtx2(a2);
		Mat3 mtx3 = mtx1 * mtx2;
		mtx3.Export(a2);
	}
	else
	{
		Mat3x4 mtx1(a1);
		Mat3x4 mtx2(a2);
		Mat3x4 mtx3 = mtx1 * mtx2;
		mtx3.Export(a2);
	}
}

double* StressStrainCppSolver::GetRadiusVector(size_t side)
{
	return _radiusVectors + side * vecStride;
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

void StressStrainCppSolver::CalculateStrainsAVX
	(
		size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		double *shiftStrains,		// ����� ����������
		double *velocityStrains,	// ����� ���. ���������
		int nodeId1,				// ����� ���� 1
		int nodeId2					// ����� ���� 2
	)
{

	// Start AVX code
	double* pmatA01 = GetRotationMatrix(nodeId1);
	double* pmatA02 = GetRotationMatrix(nodeId2);

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
	// ������� A_{21} ������������

	__m256d ivecC1 = _mm256_load_pd(GetRadiusVector(side));
	__m256d vecDP = _mm256_sub_pd(
		_mm256_load_pd(GetElementShift(nodeId1)),
		_mm256_load_pd(GetElementShift(nodeId2))); // P1-P2

	matA02el1 = _mm256_set1_pd(ivecC1.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(ivecC1.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(ivecC1.m256d_f64[2]);

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

	_mm256_store_pd(shiftStrains, res); // �������� SL, �������� ����������
	
	// ������ VL
	// ��������� ������ ������� �������� ��������� ����� ����� �2-�1 � ��1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	__declspec(align(32)) double cp1[4] = { 0 };	// ��������� ������������ 
	__declspec(align(32)) double cp2[4] = { 0 };	// ��������� ������������

	CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
	CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]
	
	__m256d cp1r = _mm256_load_pd(&cp1[0]);
	__m256d vecDV = _mm256_sub_pd(
		_mm256_load_pd(GetElementVelocity(nodeId1)),
		_mm256_load_pd(GetElementVelocity(nodeId2))); // V1-V2

	matA02el1 = _mm256_set1_pd(cp2[0]);
	matA02el2 = _mm256_set1_pd(cp2[1]);
	matA02el3 = _mm256_set1_pd(cp2[2]);

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

	_mm256_store_pd(velocityStrains, res); // �������� VL, �������� ����������

	__m256d x1 = _mm256_load_pd(GetElementShiftAngular(nodeId1));
	__m256d x2 = _mm256_load_pd(GetElementShiftAngular(nodeId2));
	_mm256_store_pd(shiftStrains + vecStride, _mm256_sub_pd(x1, x2));

	x1 = _mm256_load_pd(GetElementVelocityAngular(nodeId1));
	x2 = _mm256_load_pd(GetElementVelocityAngular(nodeId2));
	_mm256_store_pd(velocityStrains + vecStride, _mm256_sub_pd(x1, x2));
}

void StressStrainCppSolver::CalculateStrains
	(
		size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		double *shiftStrains,		// ����� ����������
		double *velocityStrains,	// ����� ���. ���������
		int nodeId1,				// ����� ���� 1
		int nodeId2					// ����� ���� 2
	)
{
	MathHelpers::Mat3x4 matA01(GetRotationMatrix(nodeId1));
	MathHelpers::Mat3x4 matA02(GetRotationMatrix(nodeId2));

	//matA01 = matA01.Tr();
	//matA02 = matA02.Tr();
	//Mat3 matA12 = matA01.Tmul(matA02);
	MathHelpers::Mat3x4 matA21 = matA02.Tmul(matA01);

	Vec3 vecC1 = MakeVec3(GetRadiusVector(side));
	Vec3 vecC2 = -vecC1;

	Vec3Ref vecP1 = MakeVec3(GetElementShift(nodeId1));
	Vec3Ref vecP2 = MakeVec3(GetElementShift(nodeId2));
	Vec3Ref vecR1 = MakeVec3(GetElementShiftAngular(nodeId1));
	Vec3Ref vecR2 = MakeVec3(GetElementShiftAngular(nodeId2));
	Vec3Ref vecV1 = MakeVec3(GetElementVelocity(nodeId1));
	Vec3Ref vecV2 = MakeVec3(GetElementVelocity(nodeId2));
	Vec3Ref vecW1 = MakeVec3(GetElementVelocityAngular(nodeId1));
	Vec3Ref vecW2 = MakeVec3(GetElementVelocityAngular(nodeId2));

	// ��������� ������ ����� ����� ����� �2-�1 � ��1
	Vec3 vecT0 = 
		vecC1 
		- matA21.Tmul(vecC2) 
		- matA01.Tmul(vecP2 - vecP1);

	vecT0.Export(shiftStrains);

	// ��������� ������ ������� �������� ��������� ����� ����� �2-�1 � ��1
	Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) 
		+ vecW1.Cross(vecC1)
		- matA21.Tmul(vecW2.Cross(vecC2));
	
	VecT1.Export(velocityStrains);

	// ��� ������� �������� ������� - ������ ������� ������� ����� � ������� ���������
	(vecR1 - vecR2).Export(shiftStrains + vecStride);
	(vecW1 - vecW2).Export(velocityStrains + vecStride);
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