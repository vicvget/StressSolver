#include "StressStrainCppSolver.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"
#include "Common.h"

#include <cstring>
#include <immintrin.h>
#include <omp.h>

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;


//#define NO_INTOMSUB


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define OUTPUT_SPACING  15
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

void StressStrainCppSolver::PrintTime()
{
	std::cout << "\n-----------------------------------\n";
	double t1 = _testTimer.Print(1, "Rotations: ");
	double t2 = _testTimer.Print(2, "Forces: ");
	double t3 = _testTimer.Print(3, "Integration: ");

	std::cout << std::setw(15) << "Summ: " << t1 + t2 + t3 << std::endl;
	//_testTimer.Print(0, "Total: ");
}

StressStrainCppSolver::StressStrainCppSolver
	(
		double* params,		// параметры
		int* links,			// индексы элементов связей (пары)
		int nLinks,			// число связей
		double *nodes,		// координаты узлов
		int nElements,		// число узлов
		double gridStep,	// шаг сетки
		double timeStep,	// шаг по времени
		int numThreads,		// число потоков
		int stride
	)
	:
		StressStrainSolver(nElements, stride),
		_isStiffnessOverriden(false),
		_isFirstSolution(true),
		//_isFirstIteration(true),
		_nIteration(0),
		_poissonRatio(0.),
		//_stageRK(0),                // этап вычисления итерации РК4 (0-3)
		//_time(0),					// время
		_timeStep(timeStep),		// шаг интегрирования
		_timeStep2(timeStep*0.5),   // половина шага интегрирования
		_timeStep4(timeStep*0.25),  // четверть шага интегрирования
		_gridStep(gridStep),		// шаг сетки
		//_elasticModulus(params[0]),	// модуль упругости
		//_dampingFactor(params[1]),	// коэффициент демпфирования
		//_density(params[2]),		// плотность материала
		_stiffScale(params[3])

{

	//_numThreads = (numThreads > 0) ? numThreads : omp_get_max_threads();

	_numThreads = (numThreads > 0) ? numThreads : 1;
	_nVariables = nElements * vecStride2;

	// масштабный коэффициент
	_gridStep2 = _gridStep * _gridStep;
	_gridStep3 = _gridStep * _gridStep2;

	// масса и момент инерции для куба
	double density = params[2];
	_cellMass = density * _gridStep3; 
	_cellInertia = _cellMass * _gridStep * _gridStep / 6;


	_elasticModulus = params[0];	// модуль упругости
	double elasticModulusScaled = _elasticModulus / _stiffScale; // отмасштабированный модуль упругости

	const double poissonFactor = 0.;
	double shearModulusScaled = elasticModulusScaled / (2 * (1 + poissonFactor));

	FindStressStrainMatrix();

	// TODO: wtf?
	double dampingFactor = params[1];	// коэффициент демпфирования
	_dampingFactorLinear = 0.01 * 2 * dampingFactor * sqrt(elasticModulusScaled * _cellMass * _gridStep);
	_dampingFactorAngular = 10 * dampingFactor * sqrt(2 * elasticModulusScaled * _cellMass * _gridStep / 3) * _gridStep2;

	_elasticFactorLinear = elasticModulusScaled * _gridStep;
	_elasticFactorAngular = 2 * shearModulusScaled * _gridStep3;

	_stressScalingFactors[0] = 1.;
	_stressScalingFactors[1] = 1.;
	_stressScalingFactors[2] = 1.;
	_puassonFactor = 0.;

	
	
	std::cout << "------------------------------" << std::endl
			  << "    CPP SOLVER IS CREATED" << std::endl
			  << "------------------------------" << std::endl
		<< std::setw(OUTPUT_SPACING) << "THREADS: " << std::setw(OUTPUT_SPACING) << _numThreads << std::endl
		<< std::setw(OUTPUT_SPACING) << "VECSTRIDE: " << std::setw(OUTPUT_SPACING) << vecStride << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFF SCALE: " << std::setw(OUTPUT_SPACING) << _stiffScale << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS L: "   << std::setw(OUTPUT_SPACING) << _elasticFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS A: " << std::setw(OUTPUT_SPACING) << _elasticFactorAngular << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING L: " << std::setw(OUTPUT_SPACING) << _dampingFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING A: " << std::setw(OUTPUT_SPACING) << _dampingFactorAngular << std::endl
		<< std::setw(OUTPUT_SPACING) << "MASS: " << std::setw(OUTPUT_SPACING) << _cellMass << std::endl
		<< std::setw(OUTPUT_SPACING) << "INERTIA: " << std::setw(OUTPUT_SPACING) << _cellInertia << std::endl
		<< std::setw(OUTPUT_SPACING) << "ELEMENTS: " << std::setw(OUTPUT_SPACING) << nElements << std::endl
		<< std::setw(OUTPUT_SPACING) << "GRID STEP: " << std::setw(OUTPUT_SPACING) << _gridStep << std::endl;

	
	// выделение памяти

	size_t varSize = _nVariables * sizeof(double);


	_initX = (double*)aligned_alloc(varSize,ALIGNMENT);
	_initDX= (double*)aligned_alloc(varSize,ALIGNMENT);
	_hDDX1 = (double*)aligned_alloc(varSize,ALIGNMENT);
	_hDDX2 = (double*)aligned_alloc(varSize,ALIGNMENT);
	_hDDX3 = (double*)aligned_alloc(varSize,ALIGNMENT);



	_varX = _dataInternal;
	_varDX = _dataInternal+_nVariables;
	_varDDX = _dataInternal+_nVariables*2;

	//_varDDX = new double[_nVariables];
	//RZ = new double[_nVariables];
	//R1Z = new double[_nVariables];
	_rotationSolver = new RotationSolver(
		nElements, 
		stride, 
		_timeStep, 
		GetElementVelocityAngular(0),
		_dataRotationMtx
		);

	// Связанные элементы: порядок x-,y-,z-,x+,y+,z+ 
	// Нумерация с 1, если 0, то связь отсутствует
	//    1
	//    |
	// 0 -x- 3
	//    |
	//    4
	// [2] = z-, [5] = z+
	_linkedElements = new int[nElements * 6];
	_radiusVectors = (double*)aligned_alloc(vecStride * 6 * sizeof(double), ALIGNMENT);
	memset(_radiusVectors, 0, vecStride * 6 * sizeof(double));

	for (int i = 0; i < 3; i++)
	{
		_radiusVectors[vecStride*3 + i*vecStride + i] = 0.5 * _gridStep;
		_radiusVectors[i*vecStride + i] = -0.5 * _gridStep;
	}

	// копируем координаты узлов сетки
	_coordinates = new double[nElements*3];
	memcpy(_coordinates, nodes, sizeof(double) * nElements * 3);
	memset(_linkedElements, 0, nElements * 6 * sizeof(int));
	for (int i = 0; i < nLinks * 2; i += 2)
	{
		// номера элементов в связях
		int elementId1 = links[i];
		int elementId2 = links[i + 1];
		int* linkedElementsOffset1 = _linkedElements + elementId1 * 6;
		int* linkedElementsOffset2 = _linkedElements + elementId2 * 6;

		// выбираем одну из координат
		double* node1 = &_coordinates[3 * links[i]];
		double* node2 = &_coordinates[3 * links[i + 1]];

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

	for (int i = 0; i < nElements; i++)
	{
		for (int dof = 0; dof < 3; dof++)
		{
			 int tmp = (_linkedElements[6 * i + dof] > 0 ? 1 : 0) + (_linkedElements[6 * i + dof + 3] > 0 ? 1 : 0);
			 _elementStressFactorCache[vecStride * i + dof] = tmp > 0 ? 1. / tmp : 0.;
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
	_testTimer.Allocate(10);
}

// virtual
StressStrainCppSolver::~StressStrainCppSolver()
{
	aligned_free(_initX);
	aligned_free(_initDX);
	aligned_free(_hDDX1);
	aligned_free(_hDDX2);
	aligned_free(_hDDX3);
	aligned_free(_radiusVectors);

	delete [] _linkedElements;
	delete [] _coordinates;
	delete _rotationSolver;
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
void StressStrainCppSolver::UpdateBuffer()
{
	GetScalarParameter(_data);
	GetVectorParameter(_dataVector);
}

float StressStrainCppSolver::UpdateBufferWithOutput()
{
	UpdateBuffer();
	float maxScalar = _data[0];
	float minScalar = _data[0];
	float avgScalar = _data[0] / _nElements;
	for (int i = 1; i < _nElements; i++)
	{
		avgScalar += _data[i] / _nElements;
		if (_data[i] > maxScalar) maxScalar = _data[i];
		if (_data[i] < minScalar) minScalar = _data[i];
	}
	std::cout << "max = " << maxScalar << " ";
	std::cout << "min = " << minScalar << " ";
	std::cout << "avg = " << avgScalar << " ";
	return maxScalar;
}

double* StressStrainCppSolver::GetRadiusVector(size_t side) const
{
	return _radiusVectors + side * vecStride;
}

int StressStrainCppSolver::GetLinkedElement(size_t elementId, size_t dof) const
{
	return _linkedElements[6 * elementId + dof];
}

void StressStrainCppSolver::OverrideScalingFactors(double stressScalingFactorX, double stressScalingFactorY, double stressScalingFactorZ)
{
	_stressScalingFactors[0] = stressScalingFactorX;
	_stressScalingFactors[1] = stressScalingFactorY;
	_stressScalingFactors[2] = stressScalingFactorZ;
}

void StressStrainCppSolver::OverrideStiffness(double elasticFactorLinear, double elasticFactorAngular, double dampingFactorLinear, double dampingFactorAngular, double stiffScale)
{
	_stiffScale = stiffScale;
	_dampingFactorAngular = dampingFactorAngular;
	_dampingFactorLinear = dampingFactorLinear;
	_elasticFactorLinear = elasticFactorLinear;
	_elasticFactorAngular = elasticFactorAngular;

	_isStiffnessOverriden = true;

	std::cout << "STIFFNESS PARAMS OVERRIDEN:" << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFF SCALE: " << std::setw(OUTPUT_SPACING) << _stiffScale << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS L: " << std::setw(OUTPUT_SPACING) << _elasticFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS A: " << std::setw(OUTPUT_SPACING) << _elasticFactorAngular << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING L: " << std::setw(OUTPUT_SPACING) << _dampingFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING A: " << std::setw(OUTPUT_SPACING) << _dampingFactorAngular << std::endl;

}

void StressStrainCppSolver::OverrideInertia(
	double mass,
	double inertia)
	{
		_cellMass = mass;
		_cellInertia = inertia;

		_isInertiaOverriden = true;

		std::cout << "INERTIA PARAMS OVERRIDEN:" << std::endl
			<< std::setw(OUTPUT_SPACING) << "MASS: " << std::setw(OUTPUT_SPACING) << _cellMass << std::endl
			<< std::setw(OUTPUT_SPACING) << "INERTIA: " << std::setw(OUTPUT_SPACING) << _cellInertia << std::endl;
	}

void StressStrainCppSolver::CrossProduct
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

//void StressStrainCppSolver::CalculateStrainsSSE
//(
//	size_t side,				// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
//	double *shiftStrains,		// выход деформаций
//	double *velocityStrains,	// выход изм. скоростей
//	size_t nodeId1,				// номер узла 1
//	size_t nodeId2				// номер узла 2
//) const
//{
//	// Start SSE code
//	double* pmatA01 = GetRotationMatrix(nodeId1);
//	double* pmatA02 = GetRotationMatrix(nodeId2);
//
//	// vecStride must be 4
//	__m128d matA01row11 = _mm_load_pd(pmatA01);
//	__m128d matA01row12 = _mm_load_pd(pmatA01 + 2);
//	__m128d matA01row21 = _mm_load_pd(pmatA01 + vecStride);
//	__m128d matA01row22 = _mm_load_pd(pmatA01 + vecStride + 2);
//	__m128d matA01row31 = _mm_load_pd(pmatA01 + vecStride2);
//	__m128d matA01row32 = _mm_load_pd(pmatA01 + vecStride2 + 2);
//
//	__m128d matA02el11;
//	__m128d matA02el12;
//	__m128d matA02el21;
//	__m128d matA02el22;
//	__m128d matA02el31;
//	__m128d matA02el32;
//
//	__m128d matA21row11;
//	__m128d matA21row12;
//	__m128d matA21row21;
//	__m128d matA21row22;
//	__m128d matA21row31;
//	__m128d matA21row32;
//
//	// matA02 column 1/3
//	matA02el11 = _mm_set1_pd(pmatA02[0]);
//
//	matA02el21 = _mm_set1_pd(pmatA02[vecStride]);
//	matA02el31 = _mm_set1_pd(pmatA02[vecStride2]);
//
//	// matA01.TMul(matA02)
//	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);
//
//	matA21row11 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
//	matA21row12 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));
//
//	// matA02 column 2/3
//	matA02el11 = _mm_set1_pd(pmatA02[1]);
//	matA02el21 = _mm_set1_pd(pmatA02[vecStride + 1]);
//	matA02el31 = _mm_set1_pd(pmatA02[vecStride2 + 1]);
//
//	// matA01.TMul(matA02)
//	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);
//
//	matA21row21 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
//	matA21row22 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));
//
//	// matA02 column 3/3
//	matA02el11 = _mm_set1_pd(pmatA02[2]);
//	matA02el21 = _mm_set1_pd(pmatA02[vecStride + 2]);
//	matA02el31 = _mm_set1_pd(pmatA02[vecStride2 + 2]);
//
//	// matA01.TMul(matA02)
//	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);
//
//	matA21row31 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
//	matA21row32 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));
//	// Матрица A_{21} сформирована
//
//	__m128d ivecC11 = _mm_load_pd(GetRadiusVector(side));
//	__m128d ivecC12 = _mm_load_pd(GetRadiusVector(side) + 2);
//	__m128d vecDP1 = _mm_sub_pd(
//		_mm_load_pd(GetElementShift(nodeId1)),
//		_mm_load_pd(GetElementShift(nodeId2))); // P1-P2
//	__m128d vecDP2 = _mm_sub_pd(
//		_mm_load_pd(GetElementShift(nodeId1) + 2),
//		_mm_load_pd(GetElementShift(nodeId2) + 2)); // P1-P2
//
//	matA02el11 = _mm_set1_pd(GetRadiusVector(side)[0]);
//	matA02el21 = _mm_set1_pd(GetRadiusVector(side)[1]);
//	matA02el31 = _mm_set1_pd(GetRadiusVector(side)[2]);
//
//	matA02el12 = _mm_mul_pd(matA21row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA21row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA21row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA21row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA21row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA21row31, matA02el31);
//
//	__m128d mul11 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
//	__m128d mul12 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);
//
//	matA02el11 = _mm_set1_pd(GetElementShift(nodeId1)[0] - GetElementShift(nodeId2)[0]);
//	matA02el21 = _mm_set1_pd(GetElementShift(nodeId1)[1] - GetElementShift(nodeId2)[1]);
//	matA02el31 = _mm_set1_pd(GetElementShift(nodeId1)[2] - GetElementShift(nodeId2)[2]);
//
//	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);
//
//	__m128d mul21 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
//	__m128d mul22 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);
//
//	__m128d res1 = _mm_add_pd(_mm_add_pd(ivecC11, mul11), mul21);
//	__m128d res2 = _mm_add_pd(_mm_add_pd(ivecC12, mul12), mul22);
//
//	_mm_store_pd(shiftStrains, res1); // получено SL, линейные компоненты
//	_mm_store_pd(shiftStrains + 2, res2);
//
//	// Расчет VL
//	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
//	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
//	// vecC2 = -vecC1
//	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)
//
//	//__declspec(align(16)) 
//	double cp1[4] = { 0 };	// векторное произведение 
//	//__declspec(align(16)) 
//	double cp2[4] = { 0 };	// векторное произведение
//
//	CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
//	CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]
//
//	__m128d cp1r1 = _mm_load_pd(&cp1[0]);
//	__m128d cp1r2 = _mm_load_pd(&cp1[2]);
//	__m128d vecDV1 = _mm_sub_pd(
//		_mm_load_pd(GetElementVelocity(nodeId1)),
//		_mm_load_pd(GetElementVelocity(nodeId2))); // V1-V2
//	__m128d vecDV2 = _mm_sub_pd(
//		_mm_load_pd(GetElementVelocity(nodeId1) + 2),
//		_mm_load_pd(GetElementVelocity(nodeId2) + 2)); // V1-V2
//
//	matA02el11 = _mm_set1_pd(cp2[0]);
//	matA02el21 = _mm_set1_pd(cp2[1]);
//	matA02el31 = _mm_set1_pd(cp2[2]);
//
//	matA02el12 = _mm_mul_pd(matA21row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA21row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA21row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA21row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA21row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA21row31, matA02el31);
//
//	mul11 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
//	mul12 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);
//
//	matA02el11 = _mm_set1_pd(GetElementVelocity(nodeId1)[0] - GetElementVelocity(nodeId2)[0]);
//	matA02el21 = _mm_set1_pd(GetElementVelocity(nodeId1)[1] - GetElementVelocity(nodeId2)[1]);
//	matA02el31 = _mm_set1_pd(GetElementVelocity(nodeId1)[2] - GetElementVelocity(nodeId2)[2]);
//
//	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
//	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
//	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
//	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
//	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
//	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);
//
//	mul21 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
//	mul22 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);
//
//	res1 = _mm_add_pd(_mm_add_pd(cp1r1, mul11), mul21);
//	res2 = _mm_add_pd(_mm_add_pd(cp1r2, mul12), mul22);
//
//	_mm_store_pd(velocityStrains, res1); // получено VL, линейные компоненты
//	_mm_store_pd(velocityStrains + 2, res2);
//
//	__m128d x11 = _mm_load_pd(GetElementShiftAngular(nodeId1));
//	__m128d x21 = _mm_load_pd(GetElementShiftAngular(nodeId2));
//	__m128d x12 = _mm_load_pd(GetElementShiftAngular(nodeId1) + 2);
//	__m128d x22 = _mm_load_pd(GetElementShiftAngular(nodeId2) + 2);
//	_mm_store_pd(shiftStrains + vecStride, _mm_sub_pd(x11, x21)); // получено VL, линейные компоненты
//	_mm_store_pd(shiftStrains + vecStride + 2, _mm_sub_pd(x12, x22));
//
//	x11 = _mm_load_pd(GetElementVelocityAngular(nodeId1));
//	x21 = _mm_load_pd(GetElementVelocityAngular(nodeId2));
//	x12 = _mm_load_pd(GetElementVelocityAngular(nodeId1) + 2);
//	x22 = _mm_load_pd(GetElementVelocityAngular(nodeId2) + 2);
//	_mm_store_pd(velocityStrains + vecStride, _mm_sub_pd(x11, x21)); // получено VL, линейные компоненты
//	_mm_store_pd(velocityStrains + vecStride + 2, _mm_sub_pd(x12, x22));
//}

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

#ifndef OMP_SOLVE
void StressStrainCppSolver::CalculateForces()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	__declspec(align(64)) double strains[8], velocityStrains[8];

	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);
	}

	const int exclusive_dofs[][2] = { { 1, 2 }, { 0, 2 }, { 1, 3 } };

	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		// обход x-,y-,z-
		for (int dof = 0; dof < 3; dof++)
		{
			size_t elementId2 = GetLinkedElement(elementId1, dof);

			if (elementId2)
			{
				elementId2--;
				CalculateStrains(dof, strains, velocityStrains, elementId1, elementId2);

				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);

				// нормальные напряжения
				GetElementStress(elementId1)[dof] += linear_strains[dof] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];
				GetElementStress(elementId2)[dof] += linear_strains[dof] * GetElementStressFactors(elementId2)[dof] * _stressScalingFactors[dof];

				// степени свободы смещений, участвующих в создании касательных напряжений
				int dof0 = exclusive_dofs[dof][0];
				int dof1 = exclusive_dofs[dof][1];

				GetElementStressAngular(elementId1)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];
				GetElementStressAngular(elementId1)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];

				GetElementStressAngular(elementId2)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId2)[dof] * _stressScalingFactors[dof];
				GetElementStressAngular(elementId2)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId2)[dof] * _stressScalingFactors[dof];

				// сила и момент из полученных деформаций
				Vec3 vForce1 = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
				Vec3 vTorque = -angular_vstrains * _dampingFactorAngular - angular_strains * _elasticFactorAngular;

				Vec3 vForce0; // сила в СК0
				Vec3 vForce2; // сила в СК второго тела

				if (vecStride == 4)
				{
					Mat3x4 matA01(GetRotationMatrix(elementId1));
					Mat3x4 matA02(GetRotationMatrix(elementId2));
					vForce0 = matA01*vForce1;
					vForce2 = matA02.Tmul(vForce0);
				}
				else
				{
					Mat3 matA01(GetRotationMatrix(elementId1));
					Mat3 matA02(GetRotationMatrix(elementId2));
					vForce0 = matA01*vForce1;
					vForce2 = matA02.Tmul(vForce0);
				}

				Vec3Ref vR = MakeVec3(GetRadiusVector(dof));
				Vec3 vForce1Torque = vR.Cross(vForce1);
				Vec3 vForce2Torque = vR.Cross(vForce2); //(-R and -vForce2 gives +vForce2Torque)

				// Full torque
				Vec3 vTorque1 = vForce1Torque + vTorque;
				Vec3 vTorque2 = vForce2Torque - vTorque;

				MakeVec3(GetElementAcceleration(elementId1)) += vForce0;
				MakeVec3(GetElementAccelerationAngular(elementId1)) += vTorque1;
				MakeVec3(GetElementAcceleration(elementId2)) -= vForce0;
				MakeVec3(GetElementAccelerationAngular(elementId2)) += vTorque2;
			}
		}
	}
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
}

#else
void StressStrainCppSolver::CalculateForces()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	__declspec(align(64)) double strains[8], velocityStrains[8];

	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);
	}

#ifdef OMP_SOLVE
#pragma omp parallel for private (strains, velocityStrains) num_threads(_numThreads)
#endif
	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{

		//memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		//memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);
		const int exclusive_dofs[][2] = { { 1, 2 }, { 0, 2 }, { 1, 3 } };

		// обход 6 связанных элементов x-,y-,z-,x+,y+,z+
		for (size_t side = 0; side < 6; side++)
		{
			int dof = side % 3;
			size_t elementId2 = _linkedElements[6 * elementId1 + side];
			if (elementId2)
			{
				elementId2--;
				CalculateStrains(side, strains, velocityStrains, elementId1, elementId2);

				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);
				size_t component = dof;
				int signFactor = (side < 3) ? 1 : -1;
				// нормальные напряжения
				GetElementStress(elementId1)[dof] += signFactor * linear_strains[dof] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];

				// степени свободы смещений, участвующих в создании касательных напряжений
				int dof0 = exclusive_dofs[dof][0];
				int dof1 = exclusive_dofs[dof][1];

				GetElementStressAngular(elementId1)[dof0] += signFactor * linear_strains[dof1] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];
				GetElementStressAngular(elementId1)[dof1] += signFactor * linear_strains[dof0] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];

				// сила и момент из полученных деформаций
				Vec3 force = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
				Vec3 torque = -angular_vstrains * _dampingFactorAngular - angular_strains * _elasticFactorAngular;

				Vec3 vAcc;
//				if (vecStride == 4)
				{
					Mat3x4 matA01(GetRotationMatrix(elementId1));
					vAcc = matA01*force;
				}
				//else
				//{
				//	Mat3 matA01(GetRotationMatrix(elementId1));
				//	vAcc = matA01*force;
				//}

				Vec3Ref vR = MakeVec3(GetRadiusVector(side));
				Vec3 forceTorque = vR.Cross(force);
				Vec3 vM = forceTorque + torque;

				MakeVec3(GetElementAcceleration(elementId1)) += vAcc;
				MakeVec3(GetElementAccelerationAngular(elementId1)) += vM;
			}
		}
	}
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
}
#endif

void StressStrainCppSolver::ApplyBoundary()
{
	//работает с типами из frm_provider
	//BCT_stresstrainBoundaryForce = 3,
	//BCT_stresstrainBoundarySealing = 4,

	double* accelerationsPointer = GetElementAcceleration(0);
	vector<BoundaryParams>::iterator it = _boundaryParamsSet.begin();

	while (it != _boundaryParamsSet.end())
	{
		switch (it->GetKind())
		{
		case 3:
			it->ApplyForceBoundary(accelerationsPointer);
			break;

		case 4:
			it->ApplySealedBoundary(accelerationsPointer);
			break;

		default:
			break;
		}
		it++;
	}
}

void StressStrainCppSolver::ApplyMass()
{
#ifdef _DEBUG
	//_controlfp(0, EM_ZERODIVIDE);
	//_control87(~_EM_ZERODIVIDE, _MCW_EM);
#endif

	//double* accelerations = GetElementAcceleration(0); // debug
	//std::cout << "Num threads = " << _numThreads << std::endl;
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		MakeVec3(GetElementAcceleration(elementId)) /= _cellMass;
		MakeVec3(GetElementAccelerationAngular(elementId)) /= _cellInertia;
	}
}

}
