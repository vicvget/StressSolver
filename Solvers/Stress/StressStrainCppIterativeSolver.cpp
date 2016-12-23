#include "StressStrainCppIterativeSolver.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

#include <cstring>

//#define NOTIMER
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define SQR(x) ((x) * (x))

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;

namespace Stress
{

StressStrainCppIterativeSolver::StressStrainCppIterativeSolver
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
		StressStrainCppSolver
			(
				params,
				links, 
				nLinks, 
				gridElements, 
				nElements, 
				gridStep,
				timeStep,
				numThreads,
				stride
			),
		_iterationNumber(0)
{
	std::cout << "NORMAL SOLVER" << std::endl << std::flush;
}

// virtual
StressStrainCppIterativeSolver::~StressStrainCppIterativeSolver()
{

}


void StressStrainCppIterativeSolver::Solve(const int nIterations)
{
	SolveFull(nIterations);
		return;
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	_iterationNumber = 0;
	_testTimer.Start(0);
	while (_iterationNumber != nIterations && _rotationSolver->IsValid())
	{
		_iterationNumber++;
		_nIteration++;
		//Solve1() used only in initialization
		Solve2();
		Solve3();
		Solve4();
		Solve5();
	}
}

void StressStrainCppIterativeSolver::SolveFull(const int nIterations)
{
//	Solve(nIterations);

	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	_iterationNumber = 0;
	_testTimer.Start(0);
	while (_iterationNumber != nIterations && _rotationSolver->IsValid())
	{
		_iterationNumber++;
		_nIteration++;

		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		_testTimer.Start(3);
		// RK4 step 1
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX1[j] = _varDDX[j] * _timeStep;
			_varX[j] += _varDX[j] * _timeStep * 0.5;
			_varDX[j] += _hDDX1[j] * 0.5;
		}
		_testTimer.Stop(3);
		//	_stageRK = 2;

		MeasuredRun(1, _rotationSolver->Solve2());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
		// RK4 step 2
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX2[j] = _varDDX[j] * _timeStep;
			_varX[j] += _hDDX1[j] * _timeStep * 0.25;
			_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
		}
		_testTimer.Stop(3);
		//	_stageRK = 3;
		MeasuredRun(1, _rotationSolver->Solve3());
		MeasuredRun(2, CalculateForces());

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		//	_time = _timeTmp + _timeStep;

		_testTimer.Start(3);
		// RK4 step 3
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			_hDDX3[j] = _varDDX[j] * _timeStep;
			_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			_varDX[j] = _initDX[j] + _hDDX3[j];
		}
		_testTimer.Stop(3);

		//	_stageRK = 4;
		MeasuredRun(1, _rotationSolver->Solve4());
		MeasuredRun(2, CalculateForces());

		// RK4 step 4
		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j++)
		{
			double sDDX = _hDDX2[j] + _hDDX3[j];
			_varX[j] = _initX[j] + (_initDX[j] + (_hDDX1[j] + sDDX) / 6.0) * _timeStep;
			_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
		}
		_testTimer.Stop(3);
		MeasuredRun(1, _rotationSolver->Solve1());
		MeasuredRun(2, CalculateForces());

		CheckVelocitySumm();
	}
#ifndef NOTIMER
	const int width = 16;
	//	_testTimer.SetWidth(width);
	std::cout << "-----------------------------------\n";
	double t1 = _testTimer.Print(1, "Rotations: ");
	double t2 = _testTimer.Print(2, "Forces: ");
	double t3 = _testTimer.Print(3, "Integration: ");

	std::cout << std::setw(width) << "Summ: " << t1 + t2 + t3 << std::endl;
	_testTimer.Print(0, "Total: ");
#endif
}

	// virtual
void StressStrainCppIterativeSolver::InitialSolve()
{
	_rotationSolver->InitialSolve();
	CalculateForces();
}


/**
* Расчет первой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_rotationSolver->Solve1();
	CalculateForces();
}

/**
* Расчет второй стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve2()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	memcpy(_initX, _varX, sizeof(double)*_nVariables);
	memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

	// RK4 step 1
	#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{
		_hDDX1[j] = _varDDX[j] * _timeStep;
		_varX[j] += _varDX[j] * _timeStep * 0.5;
		_varDX[j] += _hDDX1[j] * 0.5;       
	}
//	_stageRK = 2;
	_rotationSolver->Solve2();
	CalculateForces();
}

/**
* Расчет третьей стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve3()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	// RK4 step 2
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{      
		_hDDX2[j] = _varDDX[j] * _timeStep;
		_varX[j] += _hDDX1[j] * _timeStep * 0.25;
		_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
	}
//	_stageRK = 3;
	_rotationSolver->Solve3();
	CalculateForces();
}

/**
* Расчет четвертой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve4()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
//	_time = _timeTmp + _timeStep;

	// RK4 step 3
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{      
		_hDDX3[j] = _varDDX[j] * _timeStep;
		_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
		_varDX[j] = _initDX[j] + _hDDX3[j];
	}
//	_stageRK = 4;
	_rotationSolver->Solve4();
	CalculateForces();
}

/**
* Расчет пятой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolver::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	// RK4 step 4
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j++)
	{  
		double sDDX = _hDDX2[j]+_hDDX3[j];
		_varX[j] = _initX[j] + (_initDX[j] + (_hDDX1[j] + sDDX) / 6.0) * _timeStep;
		_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
	}
	_rotationSolver->Solve1();
	CalculateForces();

	CheckVelocitySumm();
}

/** Получить смещения
* @param data - массив для записи смещений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetDisplacement
	(
		float* data
	)
{
	int internalIndex = 0;
	double displacement;
	double squareSum;

	for(size_t elementId = 0; elementId < _nElements; elementId++)
	{		
		squareSum = 0;
		for (int j = 0; j < 3; j++)
		{
			displacement = GetElementGridCoordinates(elementId)[j] - GetElementShift(elementId)[j];
			squareSum += SQR(displacement);
		}
		data[elementId] = (float)sqrt(squareSum);
	}
}

/** Получить напряжения по первой теории прочности
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesByFirstTheoryOfStrength
	(
		float* data
	)
{
	double* movements = _dataInternal;
	double coordinates[3];
	double relativeShifts[3];
	double relativeShiftsSigned[3];
	int neighbourNodeNumber;
	int shiftIndex;
	double shift;
	double fullRelativeShift;
	double maxRelativeShift = 0;

	for (int i = 0; i < _nElements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			coordinates[j] = _dataInternal[6 * i + j];
			relativeShifts[j] = 0;
			relativeShiftsSigned[j] = 0;
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				neighbourNodeNumber = _linkedElements[6 * i + j + 3 * k];
				if (neighbourNodeNumber != 0)
				{
					neighbourNodeNumber--;
					shiftIndex = j;
					shift = movements[6 * neighbourNodeNumber + shiftIndex] - coordinates[shiftIndex];
					if (fabs(shift) > fabs(relativeShifts[shiftIndex]))
					{
						relativeShiftsSigned[shiftIndex] = shift * (2 * k - 1) - _gridStep;
						if (fabs(relativeShiftsSigned[shiftIndex]) < 1e-15)
						{
							relativeShiftsSigned[shiftIndex] = 0.0;
						}
						relativeShifts[shiftIndex] = fabs(shift);
					}
				}
			}
		}
		fullRelativeShift  = relativeShiftsSigned[0];
		for (int j = 1; j < 3; j++)
		{
			if (fabs(relativeShiftsSigned[j]) > fabs(fullRelativeShift))
			{
				fullRelativeShift = relativeShiftsSigned[j];
			}
		}
		data[i] = (float)fabs(_elasticModulus * fullRelativeShift / _gridStep);
	}
}


/** Получить напряжения по X
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesX(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		data[elementId] = (float)(GetElementStress(elementId)[0] * _elasticFactorLinear / (_gridStep *_gridStep));
	}
}

/** Получить напряжения по Y
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesY(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		data[elementId] = (float)(GetElementStress(elementId)[1] * _elasticFactorLinear / (_gridStep *_gridStep));
	}
}

/** Получить напряжения по Z
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesZ(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		data[elementId] = (float)(GetElementStress(elementId)[2] * _elasticFactorLinear / (_gridStep *_gridStep));
	}
}

/** Получить напряжения по XY
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesXY(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		// TODO: scalingFactorX и scalingFactorY для 2 компонент!
		data[elementId] = (float)(GetElementStressAngular(elementId)[2] * _elasticFactorLinear / 2. /(1.+0.28) / (_gridStep * _gridStep));
	}
}

/** Получить напряжения по XZ
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesXZ(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		data[elementId] = (float)(GetElementStressAngular(elementId)[1] * _elasticFactorLinear / 2. / (1. + 0.28) / (_gridStep * _gridStep)); // костыль для толщины в 10 мм
	}
}

/** Получить напряжения по YZ
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesYZ(float* data)
{
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		data[elementId] = (float)(GetElementStressAngular(elementId)[0] * _elasticFactorLinear / 2. / (1. + 0.28) / (_gridStep * _gridStep)); // костыль для толщины в 10 мм
	}
}



/** Получить напряжения по von Mises
* @param data - массив для записи напряжений как скалярного параметра
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesByVonMises
	(
		float* data
	)
{
	float* sigma1 = data + _nElements;
	float* sigma2 = data + _nElements*2;
	float* sigma3 = data + _nElements*3;
	float* sigma4 = data + _nElements*4;
	float* sigma5 = data + _nElements*5;
	float* sigma6 = data + _nElements*6;

#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{	
		data[elementId] = (float)sqrt
			(
					0.5 *
					(SQR(sigma1[elementId] - sigma2[elementId]) +
					SQR(sigma2[elementId] - sigma3[elementId]) +
					SQR(sigma1[elementId] - sigma3[elementId]) +
					6 * 
					(SQR(sigma4[elementId]) +
					SQR(sigma5[elementId]) +
					SQR(sigma6[elementId])
					))
			);
	}
}

void StressStrainCppIterativeSolver::CalculateStrains
(
size_t side,				// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *shiftStrains,		// выход деформаций
double *velocityStrains,	// выход изм. скоростей
size_t nodeId1,				// номер узла 1
size_t nodeId2				// номер узла 2
) const
{
	//CalculateStrainsAVX(side, shiftStrains, velocityStrains, nodeId1, nodeId2);
	//return;
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

void StressStrainCppIterativeSolver::CalculateStrainsUa
(
size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
double *shiftStrains,		// выход деформаций
double *velocityStrains,	// выход изм. скоростей
size_t nodeId1,				// номер узла 1
size_t nodeId2					// номер узла 2
) const
{
	MathHelpers::Mat3 matA01(GetRotationMatrix(nodeId1));
	MathHelpers::Mat3 matA02(GetRotationMatrix(nodeId2));

	//matA01 = matA01.Tr();
	//matA02 = matA02.Tr();
	//Mat3 matA12 = matA01.Tmul(matA02);
	MathHelpers::Mat3 matA21 = matA02.Tmul(matA01);

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

	//PrintVector(vecC1,"vecC1");
	//PrintVector(vecP1, "vecP1");
	//PrintVector(vecP2,"vecP2");
	//PrintVector(vecR1,"vecR1");
	//PrintVector(vecR2,"vecR2");
	//PrintVector(vecV1, "vecV1");
	//PrintVector(vecV2,"vecV2");
	//PrintMatrix(matA01, "matA01");
	//PrintMatrix(matA02, "matA02");
	//PrintMatrix(matA21, "matA21");

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


}