#include "StressStrainCppSolver.h"

#include "CsrSymmetricMatrix.h"
#include "../../AdditionalModules/fmath/Vector3.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"

#include <omp.h>
#include "../../AdditionalModules/fmath/Matrix3x4.h"
#include "Common.h"


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
	_numThreads = (numThreads > 0) ? numThreads : omp_get_max_threads();

	_nVariables = nElements * vecStride2;

	// масштабный коэффициент
	_gridStep2 = _gridStep * _gridStep;
	_gridStep3 = _gridStep * _gridStep2;

	// масса и момент инерции для куба
	double density = params[2];
	_cellMass = density * _gridStep3; 
	_cellInertia = _cellMass * _gridStep * _gridStep / 6;


	double _elasticModulus = params[0];	// модуль упругости
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


	// DEBUG!!
	//_shearModulusScaled = _elasticModulusScaled / 100;
	//_dampingFactorLinear = 10;
	//_dampingFactorAngular = 10;
	
	
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
	_elements = new double[nElements*3];
	memcpy(_elements, nodes, sizeof(double) * nElements * 3);
	memset(_linkedElements, 0, nElements * 6 * sizeof(int));
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
	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_dataInternal[i * 2 * vecStride + j] = nodes[i * 3 + j];
			if (i == 2 && j == 0)
				_dataInternal[i * 2 * vecStride + j] += 0.01;
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
	delete [] _elements;
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


double* StressStrainCppSolver::GetRadiusVector(size_t side) const
{
	return _radiusVectors + side * vecStride;
}

void StressStrainCppSolver::OverrideStiffness(double elasticFactorLinear, double elasticFactorAngular, double dampingFactorLinear, double dampingFactorAngular, double stiffScale)
{
	_stiffScale = stiffScale;
	_dampingFactorAngular = dampingFactorAngular;
	_dampingFactorLinear = dampingFactorLinear;
	_elasticFactorLinear = elasticFactorLinear;
	_elasticFactorAngular = elasticFactorAngular;

	_isStiffnessOverriden = true;

	std::cout << "PARAMS OVERRIDEN:" << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFF SCALE: " << std::setw(OUTPUT_SPACING) << _stiffScale << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS L: " << std::setw(OUTPUT_SPACING) << _elasticFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "STIFFNESS A: " << std::setw(OUTPUT_SPACING) << _elasticFactorAngular << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING L: " << std::setw(OUTPUT_SPACING) << _dampingFactorLinear << std::endl
		<< std::setw(OUTPUT_SPACING) << "DAMPING A: " << std::setw(OUTPUT_SPACING) << _dampingFactorAngular << std::endl;

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


void StressStrainCppSolver::CalculateStrainsSSE
(
	size_t side,				// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
	double *shiftStrains,		// выход деформаций
	double *velocityStrains,	// выход изм. скоростей
	size_t nodeId1,				// номер узла 1
	size_t nodeId2				// номер узла 2
) const
{
	// Start SSE code
	double* pmatA01 = GetRotationMatrix(nodeId1);
	double* pmatA02 = GetRotationMatrix(nodeId2);

	// vecStride must be 4
	__m128d matA01row11 = _mm_load_pd(pmatA01);
	__m128d matA01row12 = _mm_load_pd(pmatA01 + 2);
	__m128d matA01row21 = _mm_load_pd(pmatA01 + vecStride);
	__m128d matA01row22 = _mm_load_pd(pmatA01 + vecStride + 2);
	__m128d matA01row31 = _mm_load_pd(pmatA01 + vecStride2);
	__m128d matA01row32 = _mm_load_pd(pmatA01 + vecStride2 + 2);

	__m128d matA02el11;
	__m128d matA02el12;
	__m128d matA02el21;
	__m128d matA02el22;
	__m128d matA02el31;
	__m128d matA02el32;

	__m128d matA21row11;
	__m128d matA21row12;
	__m128d matA21row21;
	__m128d matA21row22;
	__m128d matA21row31;
	__m128d matA21row32;

	// matA02 column 1/3
	matA02el11 = _mm_set1_pd(pmatA02[0]);

	matA02el21 = _mm_set1_pd(pmatA02[vecStride]);
	matA02el31 = _mm_set1_pd(pmatA02[vecStride2]);

	// matA01.TMul(matA02)
	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);

	matA21row11 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
	matA21row12 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));

	// matA02 column 2/3
	matA02el11 = _mm_set1_pd(pmatA02[1]);
	matA02el21 = _mm_set1_pd(pmatA02[vecStride + 1]);
	matA02el31 = _mm_set1_pd(pmatA02[vecStride2 + 1]);

	// matA01.TMul(matA02)
	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);

	matA21row21 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
	matA21row22 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));

	// matA02 column 3/3
	matA02el11 = _mm_set1_pd(pmatA02[2]);
	matA02el21 = _mm_set1_pd(pmatA02[vecStride + 2]);
	matA02el31 = _mm_set1_pd(pmatA02[vecStride2 + 2]);

	// matA01.TMul(matA02)
	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);

	matA21row31 = _mm_add_pd(matA02el11, _mm_add_pd(matA02el21, matA02el31));
	matA21row32 = _mm_add_pd(matA02el12, _mm_add_pd(matA02el22, matA02el32));
	// Матрица A_{21} сформирована

	__m128d ivecC11 = _mm_load_pd(GetRadiusVector(side));
	__m128d ivecC12 = _mm_load_pd(GetRadiusVector(side) + 2);
	__m128d vecDP1 = _mm_sub_pd(
		_mm_load_pd(GetElementShift(nodeId1)),
		_mm_load_pd(GetElementShift(nodeId2))); // P1-P2
	__m128d vecDP2 = _mm_sub_pd(
		_mm_load_pd(GetElementShift(nodeId1) + 2),
		_mm_load_pd(GetElementShift(nodeId2) + 2)); // P1-P2

	matA02el11 = _mm_set1_pd(GetRadiusVector(side)[0]);
	matA02el21 = _mm_set1_pd(GetRadiusVector(side)[1]);
	matA02el31 = _mm_set1_pd(GetRadiusVector(side)[2]);

	matA02el12 = _mm_mul_pd(matA21row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA21row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA21row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA21row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA21row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA21row31, matA02el31);

	__m128d mul11 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
	__m128d mul12 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);

	matA02el11 = _mm_set1_pd(GetElementShift(nodeId1)[0] - GetElementShift(nodeId2)[0]);
	matA02el21 = _mm_set1_pd(GetElementShift(nodeId1)[1] - GetElementShift(nodeId2)[1]);
	matA02el31 = _mm_set1_pd(GetElementShift(nodeId1)[2] - GetElementShift(nodeId2)[2]);

	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);

	__m128d mul21 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
	__m128d mul22 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);

	__m128d res1 = _mm_add_pd(_mm_add_pd(ivecC11, mul11), mul21);
	__m128d res2 = _mm_add_pd(_mm_add_pd(ivecC12, mul12), mul22);

	_mm_store_pd(shiftStrains, res1); // получено SL, линейные компоненты
	_mm_store_pd(shiftStrains + 2, res2);

	// Расчет VL
	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	//__declspec(align(16)) 
	double cp1[4] = { 0 };	// векторное произведение 
	//__declspec(align(16)) 
	double cp2[4] = { 0 };	// векторное произведение

	CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
	CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]

	__m128d cp1r1 = _mm_load_pd(&cp1[0]);
	__m128d cp1r2 = _mm_load_pd(&cp1[2]);
	__m128d vecDV1 = _mm_sub_pd(
		_mm_load_pd(GetElementVelocity(nodeId1)),
		_mm_load_pd(GetElementVelocity(nodeId2))); // V1-V2
	__m128d vecDV2 = _mm_sub_pd(
		_mm_load_pd(GetElementVelocity(nodeId1) + 2),
		_mm_load_pd(GetElementVelocity(nodeId2) + 2)); // V1-V2

	matA02el11 = _mm_set1_pd(cp2[0]);
	matA02el21 = _mm_set1_pd(cp2[1]);
	matA02el31 = _mm_set1_pd(cp2[2]);

	matA02el12 = _mm_mul_pd(matA21row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA21row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA21row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA21row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA21row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA21row31, matA02el31);

	mul11 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
	mul12 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);

	matA02el11 = _mm_set1_pd(GetElementVelocity(nodeId1)[0] - GetElementVelocity(nodeId2)[0]);
	matA02el21 = _mm_set1_pd(GetElementVelocity(nodeId1)[1] - GetElementVelocity(nodeId2)[1]);
	matA02el31 = _mm_set1_pd(GetElementVelocity(nodeId1)[2] - GetElementVelocity(nodeId2)[2]);

	matA02el12 = _mm_mul_pd(matA01row12, matA02el11);
	matA02el11 = _mm_mul_pd(matA01row11, matA02el11);
	matA02el22 = _mm_mul_pd(matA01row22, matA02el21);
	matA02el21 = _mm_mul_pd(matA01row21, matA02el21);
	matA02el32 = _mm_mul_pd(matA01row32, matA02el31);
	matA02el31 = _mm_mul_pd(matA01row31, matA02el31);

	mul21 = _mm_add_pd(_mm_add_pd(matA02el11, matA02el21), matA02el31);
	mul22 = _mm_add_pd(_mm_add_pd(matA02el12, matA02el22), matA02el32);

	res1 = _mm_add_pd(_mm_add_pd(cp1r1, mul11), mul21);
	res2 = _mm_add_pd(_mm_add_pd(cp1r2, mul12), mul22);

	_mm_store_pd(velocityStrains, res1); // получено VL, линейные компоненты
	_mm_store_pd(velocityStrains + 2, res2);

	__m128d x11 = _mm_load_pd(GetElementShiftAngular(nodeId1));
	__m128d x21 = _mm_load_pd(GetElementShiftAngular(nodeId2));
	__m128d x12 = _mm_load_pd(GetElementShiftAngular(nodeId1) + 2);
	__m128d x22 = _mm_load_pd(GetElementShiftAngular(nodeId2) + 2);
	_mm_store_pd(shiftStrains + vecStride, _mm_sub_pd(x11, x21)); // получено VL, линейные компоненты
	_mm_store_pd(shiftStrains + vecStride + 2, _mm_sub_pd(x12, x22));

	x11 = _mm_load_pd(GetElementVelocityAngular(nodeId1));
	x21 = _mm_load_pd(GetElementVelocityAngular(nodeId2));
	x12 = _mm_load_pd(GetElementVelocityAngular(nodeId1) + 2);
	x22 = _mm_load_pd(GetElementVelocityAngular(nodeId2) + 2);
	_mm_store_pd(velocityStrains + vecStride, _mm_sub_pd(x11, x21)); // получено VL, линейные компоненты
	_mm_store_pd(velocityStrains + vecStride + 2, _mm_sub_pd(x12, x22));
}

void StressStrainCppSolver::CalculateStrainsAVX
	(
		size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		double *shiftStrains,		// выход деформаций
		double *velocityStrains,	// выход изм. скоростей
		size_t nodeId1,				// номер узла 1
		size_t nodeId2					// номер узла 2
		) const
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
	// Матрица A_{21} сформирована

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

	_mm256_store_pd(shiftStrains, res); // получено SL, линейные компоненты
	
	// Расчет VL
	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	__declspec(align(32)) double cp1[4] = { 0 };	// векторное произведение 
	__declspec(align(32)) double cp2[4] = { 0 };	// векторное произведение

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

	_mm256_store_pd(velocityStrains, res); // получено VL, линейные компоненты

	__m256d x1 = _mm256_load_pd(GetElementShiftAngular(nodeId1));
	__m256d x2 = _mm256_load_pd(GetElementShiftAngular(nodeId2));
	_mm256_store_pd(shiftStrains + vecStride, _mm256_sub_pd(x1, x2));

	x1 = _mm256_load_pd(GetElementVelocityAngular(nodeId1));
	x2 = _mm256_load_pd(GetElementVelocityAngular(nodeId2));
	_mm256_store_pd(velocityStrains + vecStride, _mm256_sub_pd(x1, x2));
}

void StressStrainCppSolver::CalculateStrainsFMA
(
size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *shiftStrains,		// выход деформаций
double *velocityStrains,	// выход изм. скоростей
size_t nodeId1,				// номер узла 1
size_t nodeId2					// номер узла 2
) const
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
	matA21row1 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row1 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row1);
	matA21row1 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row1);

	// matA02 column 2/3
	matA02el1 = _mm256_set1_pd(pmatA02[1]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 1]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 1]);

	// matA01.TMul(matA02)
	matA21row2 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row2 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row2);
	matA21row2 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row2);

	// matA02 column 3/3
	matA02el1 = _mm256_set1_pd(pmatA02[2]);
	matA02el2 = _mm256_set1_pd(pmatA02[vecStride + 2]);
	matA02el3 = _mm256_set1_pd(pmatA02[vecStride2 + 2]);

	// matA01.TMul(matA02)
	matA21row3 = _mm256_mul_pd(matA01row1, matA02el1);
	matA21row3 = _mm256_fmadd_pd(matA01row2, matA02el2, matA21row3);
	matA21row3 = _mm256_fmadd_pd(matA01row3, matA02el3, matA21row3);
	// Матрица A_{21} сформирована

	__m256d ivecC1 = _mm256_load_pd(GetRadiusVector(side));
	__m256d vecDP = _mm256_sub_pd(
		_mm256_load_pd(GetElementShift(nodeId1)),
		_mm256_load_pd(GetElementShift(nodeId2))); // P1-P2

	matA02el1 = _mm256_set1_pd(ivecC1.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(ivecC1.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(ivecC1.m256d_f64[2]);

	__m256d 
	res = _mm256_fmadd_pd(matA21row1, matA02el1, ivecC1);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);

	matA02el1 = _mm256_set1_pd(vecDP.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(vecDP.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(vecDP.m256d_f64[2]);

	res = _mm256_fmadd_pd(matA01row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA01row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA01row3, matA02el3, res);

	_mm256_store_pd(shiftStrains, res); // получено SL, линейные компоненты

	// Расчет VL
	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
	// vecC2 = -vecC1
	// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

	__declspec(align(32)) double cp1[4] = { 0 };	// векторное произведение 
	__declspec(align(32)) double cp2[4] = { 0 };	// векторное произведение

	CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
	CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]

	res = _mm256_load_pd(&cp1[0]);
	__m256d vecDV = _mm256_sub_pd(
		_mm256_load_pd(GetElementVelocity(nodeId1)),
		_mm256_load_pd(GetElementVelocity(nodeId2))); // V1-V2

	matA02el1 = _mm256_set1_pd(cp2[0]);
	matA02el2 = _mm256_set1_pd(cp2[1]);
	matA02el3 = _mm256_set1_pd(cp2[2]);

	res = _mm256_fmadd_pd(matA21row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA21row2, matA02el2, res);

	matA02el1 = _mm256_set1_pd(vecDV.m256d_f64[0]);
	matA02el2 = _mm256_set1_pd(vecDV.m256d_f64[1]);
	matA02el3 = _mm256_set1_pd(vecDV.m256d_f64[2]);

	res = _mm256_fmadd_pd(matA01row1, matA02el1, res);
	res = _mm256_fmadd_pd(matA01row2, matA02el2, res);
	res = _mm256_fmadd_pd(matA01row3, matA02el3, res);


	_mm256_store_pd(velocityStrains, res); // получено VL, линейные компоненты

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
		double *shiftStrains,		// выход деформаций
		double *velocityStrains,	// выход изм. скоростей
		size_t nodeId1,				// номер узла 1
		size_t nodeId2					// номер узла 2
		) const
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

	// переводим вектор линии точек связи С2-С1 в СК1
	Vec3 vecT0 = 
		vecC1 
		- matA21.Tmul(vecC2) 
		- matA01.Tmul(vecP2 - vecP1);

	vecT0.Export(shiftStrains);

	// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
	Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) 
		+ vecW1.Cross(vecC1)
		- matA21.Tmul(vecW2.Cross(vecC2));
	
	VecT1.Export(velocityStrains);

	// для угловых степеней свободы - просто берутся разница углов и угловых скоростей
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