#include "StressStrainCppIterativeSolver.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

#include <cstring>

#define NOTIMER
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
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	FTimer test_timer;
	//std::cout << "Full routine\n";
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

void StressStrainCppIterativeSolver::SolveFull(const int nIteratons)
{
	Solve(nIteratons);
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
		data[elementId] = GetElementStress(elementId)[0] * _elasticFactorLinear / (_gridStep *_gridStep) * _stressScalingFactorX;
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
		data[elementId] = GetElementStress(elementId)[1] * _elasticFactorLinear / (_gridStep *_gridStep) * _stressScalingFactorY;
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
		data[elementId] = GetElementStress(elementId)[2] * _elasticFactorLinear / (_gridStep *_gridStep) * _stressScalingFactorZ;
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
		data[elementId] = GetElementStressAngular(elementId)[2] * _elasticFactorLinear / 2. / (_gridStep * _gridStep) * _stressScalingFactorX;
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
		data[elementId] = GetElementStressAngular(elementId)[1] * _elasticFactorLinear / 2. / (_gridStep * 0.01); // костыль для толщины в 10 мм
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
		data[elementId] = GetElementStressAngular(elementId)[0] * _elasticFactorLinear / 2. / (_gridStep * 0.01); // костыль для толщины в 10 мм
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
	double relativeShiftsSigned[6];

	double sigma[6];
	
#pragma omp parallel for private(relativeShiftsSigned, sigma) num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		
		for (size_t dof = 0; dof < 3; dof++)
		{
			relativeShiftsSigned[dof] = GetElementStress(elementId)[dof];
			relativeShiftsSigned[dof] /= _gridStep;
			relativeShiftsSigned[dof+3] = GetElementStressAngular(elementId)[dof];
		}

		double summsq = 0;

		for (int i = 0; i < 6; i++)
		{
			sigma[i] = 0;
			for (int j = 0; j < 6; j++)
			{
				sigma[i] += _lameMatrix[i][j] * relativeShiftsSigned[j];
			}
			if (i > 2)
			{
				sigma[i] *= 2;
				summsq += SQR(sigma[i]);
			}
		}
		//summsq = 0;
		data[elementId] = (float)(_elasticModulus * sqrt
			(
				0.5 *
					(
						SQR(sigma[0] - sigma[1]) +
						SQR(sigma[1] - sigma[2]) +
						SQR(sigma[0] - sigma[2]) +
						6 * summsq
					)
			));
	}
}

void StressStrainCppIterativeSolver::CalculateForces()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	__declspec(align(32)) double strains[8], velocityStrains[8];
#if 1
	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);
	}

	const int exclusive_dofs[][2] = { { 1, 2 }, { 0, 2 }, { 1, 3 } };

	//#pragma omp parallel for private (strains, velocityStrains) num_threads(_numThreads)
	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		// обход x-,y-,z-
		double* accelerationVector1 = GetElementAcceleration(elementId1);
		double* stressVector1 = GetElementStress(elementId1);
		for (int dof = 0; dof < 3; dof++)
		{
			size_t elementId2 = GetLinkedElement(elementId1, dof);

			if (elementId2)
			{
				elementId2--;
				CalculateStrainsAVX(dof, strains, velocityStrains, elementId1, elementId2);

				double* accelerationVector2 = GetElementAcceleration(elementId1);
				double* stressVector2 = GetElementStress(elementId1);


				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);

				// нормальные напряжения
				GetElementStress(elementId1)[dof] += linear_strains[dof] * GetElementStressFactors(elementId1)[dof];
				GetElementStress(elementId2)[dof] += linear_strains[dof] * GetElementStressFactors(elementId2)[dof];

				// касательные ??? - todo: переделать на сдвиг
				int dof0 = exclusive_dofs[dof][0];
				int dof1 = exclusive_dofs[dof][1];
				GetElementStressAngular(elementId1)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId1)[dof];
				GetElementStressAngular(elementId1)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId1)[dof];
				GetElementStressAngular(elementId2)[dof0] += linear_strains[dof1] * GetElementStressFactors(elementId2)[dof];
				GetElementStressAngular(elementId2)[dof1] += linear_strains[dof0] * GetElementStressFactors(elementId2)[dof];
			
				//GetElementStressAngular(elementId1)[dof] += angular_strains[dof] * GetElementStressFactors(elementId1)[dof];
				//GetElementStressAngular(elementId2)[dof] += angular_strains[dof] * GetElementStressFactors(elementId2)[dof];

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
#else
	//#pragma omp parallel for private (strains, velocityStrains) num_threads(_numThreads)
	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{

		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);

		_testTimer.Start(5);

		// обход 6 связанных элементов x-,y-,z-,x+,y+,z+
		
		int dofFactors[3] = {0};


		for (size_t side = 0; side < 6; side++)
		{
			int dof = side % 3;
			dofFactors[dof]++;
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
				int testElement = 55;
				if (_nIteration == 500)
				{
					if (elementId1 == testElement && dof == 0)
					{
						if (side < 3)
						{
							std::cout << "1SFb=" << dofFactors[dof] << " S=" << GetElementStress(elementId1)[dof] << " Linear=" << linear_strains[dof] << std::endl;
						}
						else
						{
							std::cout << "2SFb=" << dofFactors[dof] << " S=" << GetElementStress(elementId1)[dof] << " Linear=" << linear_strains[dof] << std::endl;
						}
					}
				}

				//for (size_t component = 0; component < 3; component++)
				{
					if (side < 3)
					{
						GetElementStress(elementId1)[component] += linear_strains[component];
						GetElementStressAngular(elementId1)[component] += angular_strains[component];
					}
					else
					{
						GetElementStress(elementId1)[component] -= linear_strains[component];
						GetElementStressAngular(elementId1)[component] -= angular_strains[component];						
						GetElementStress(elementId1)[component] /= dofFactors[dof];
					}
				}

				if (_nIteration == 500)
				{
					if (elementId1 == testElement && dof == 0)
					{
						if(side < 3)
						{
							std::cout << "1SFa=" << dofFactors[dof] << " S=" << GetElementStress(elementId1)[dof] << " Linear=" << linear_strains[dof] << std::endl;
						}
						else
						{
							std::cout << "2SFa=" << dofFactors[dof] << " S=" << GetElementStress(elementId1)[dof] << " Linear=" << linear_strains[dof] << std::endl;
						}

					}
				}


				// сила и момент из полученных деформаций
				Vec3 force  = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
				Vec3 torque = -angular_vstrains * _dampingFactorAngular - angular_strains * _elasticFactorAngular;

				// отладочный вывод деформаций
				df[0] = -strains[0];
				df[1] = -strains[1];
				df[2] = -strains[2];
				df[3] = -strains[0 + vecStride];
				df[4] = -strains[1 + vecStride];
				df[5] = -strains[2 + vecStride];
				df[6] = force.X();
				df[7] = force.Y();
				df[8] = force.Z();
				df[9] = torque.X();
				df[10] = torque.Y();
				df[11] = torque.Z();


				Vec3 vAcc;
				if (vecStride == 4)
				{
					Mat3x4 matA01(GetRotationMatrix(elementId1));
					vAcc = matA01*force;
				}
				else
				{
					Mat3 matA01(GetRotationMatrix(elementId1));
					vAcc = matA01*force;
				}

				Vec3Ref vR = MakeVec3(GetRadiusVector(side));
				Vec3 forceTorque = vR.Cross(force);
				Vec3 vM = forceTorque + torque;

				MakeVec3(GetElementAcceleration(elementId1)) += vAcc;
				MakeVec3(GetElementAccelerationAngular(elementId1)) += vM;
			}
		}
		_testTimer.Stop(5);

		//Vec3Ref force = MakeVec3(&forces[0]);
		//Vec3Ref torque = MakeVec3(&forces[0] + vecStride);
		
	}
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
#endif
}


void StressStrainCppIterativeSolver::ApplyBoundary()
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

void StressStrainCppIterativeSolver::ApplyMass()
{
#ifdef _DEBUG
	//_controlfp(0, EM_ZERODIVIDE);
	//_control87(~_EM_ZERODIVIDE, _MCW_EM);
#endif

	//double* accelerations = GetElementAcceleration(0); // debug
	//std::cout << "Num threads = " << _numThreads << std::endl;
#pragma omp parallel for num_threads(_numThreads)
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		MakeVec3(GetElementAcceleration(elementId)) /= _cellMass;
		MakeVec3(GetElementAccelerationAngular(elementId)) /= _cellInertia;
	}
}
}