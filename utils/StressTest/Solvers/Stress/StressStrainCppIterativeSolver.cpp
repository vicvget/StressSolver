#include "StressStrainCppIterativeSolver.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"


//#define NOLINKSH
//#define NO_INTOMSUB
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
* ������ ������ ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolver::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_rotationSolver->Solve1();
	CalculateForces();
}

/**
* ������ ������ ������ ������ �����-�����
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
* ������ ������� ������ ������ �����-�����
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
* ������ ��������� ������ ������ �����-�����
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
* ������ ����� ������ ������ �����-�����
*/
// virtual
void StressStrainCppIterativeSolver::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	for (int j = 0; j < _nVariables; j++)
	{  
		double sDDX = _hDDX2[j]+_hDDX3[j];
		_varX[j] = _initX[j] + (_initDX[j] + (_hDDX1[j] + sDDX) / 6.0) * _timeStep;
		_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
	}
	_rotationSolver->Solve1();
	CalculateForces();
}

/** �������� ��������
* @param data - ������ ��� ������ �������� ��� ���������� ���������
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

/** �������� ���������� �� ������ ������ ���������
* @param data - ������ ��� ������ ���������� ��� ���������� ���������
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

/** �������� ���������� �� von Mises
* @param data - ������ ��� ������ ���������� ��� ���������� ���������
*/
// virtual
void StressStrainCppIterativeSolver::GetStressesByVonMises
	(
		float* data
	)
{
	double relativeShiftsSigned[6];

	double sigma[6];
	
	for (size_t elementId = 0; elementId < _nElements; elementId++)
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

	for (size_t elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		
		double* accelerationVector = GetElementAcceleration(elementId1);

		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);

		_testTimer.Start(5);

		// ����� 6 ��������� ��������� x-,y-,z-,x+,y+,z+
		for (size_t side = 0; side < 6; side++)
		{
			size_t elementId2 = _linkedElements[6 * elementId1 + side];
			if (elementId2)
			{
				elementId2--;
#ifdef ALIGNED_MEM
				CalculateStrainsAVX(side, strains, velocityStrains, elementId1, elementId2);
#endif
///				CalculateStrains(side, strains, velocityStrains, elementId1, elementId2);

				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);

				for (size_t component = 0; component < 3; component++)
				{
					GetElementStress(elementId1)[component] += strains[component];
					GetElementStress(elementId1)[component + vecStride] += strains[component + vecStride];
				}

				// ���� � ������ �� ���������� ����������
				Vec3 force  = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
				Vec3 torque = -angular_vstrains * _dampingFactorAngular - angular_strains * _elasticFactorAngular;

				// ���������� ����� ����������
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

				//MakeVec3(GetElementAcceleration(elementId1)) += vAcc;
				//MakeVec3(GetElementAccelerationAngular(elementId1)) += vM;
				
				
				for (int i = 0; i < 3; i++)
				{
					accelerationVector[i] += vAcc[i];
					accelerationVector[i + vecStride] += vM[i];
				}

			}
		}
		_testTimer.Stop(5);

		//Vec3Ref force = MakeVec3(&forces[0]);
		//Vec3Ref torque = MakeVec3(&forces[0] + vecStride);
		
	}
	ApplyBoundary(); // ������������ ���� � �������
	ApplyMass();	 // ��������� ��������� �������� ��� �� ����� � �������� �� ������� �������
}


void StressStrainCppIterativeSolver::ApplyBoundary()
{
	//�������� � ������ �� frm_provider
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

	double* accelerations = GetElementAcceleration(0); // debug
	for (int elementId = 0; elementId < _nElements; elementId++)
	{
		MakeVec3(GetElementAcceleration(elementId)) /= _cellMass;
		MakeVec3(GetElementAccelerationAngular(elementId)) /= _cellInertia;
	}
}
}