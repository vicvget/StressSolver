#ifndef USE_KNC

#include "StressStrainCppIterativeSolverAVX.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;

namespace Stress
{

StressStrainCppIterativeSolverAVX::StressStrainCppIterativeSolverAVX
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
		StressStrainCppIterativeSolver
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
			)
{
	std::cout << "AVX SOLVER" << std::endl << std::flush;
}

// virtual
StressStrainCppIterativeSolverAVX::~StressStrainCppIterativeSolverAVX()
{

}

// AVX-версия
void StressStrainCppIterativeSolverAVX::SolveFull(const int nIterations)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep2;
	__m256d timeStep4;
	__m256d constantD2;
	__m256d constantD6;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep2 = _mm256_set1_pd(_timeStep2);
	timeStep4 = _mm256_set1_pd(_timeStep4);
	constantD2 = _mm256_set1_pd(0.5);
	constantD6 = _mm256_set1_pd(1 / 6.0);


	FTimer test_timer;
	_iterationNumber = 0;
	_testTimer.Start(0);

	__m256d Xtmp, DXtmp, DDXtmp, tmp;
	__m256d hDDX1, hDDX2, hDDX3, sDDX;

	while (_iterationNumber != nIterations && _rotationSolver->IsValid())
	{
		_iterationNumber++;
		_nIteration++;

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		_testTimer.Start(3);
	#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX1[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _varDX[j] * _timeStep2;
			//_varDX[j] += _hDDX1[j] * 0.5;

			Xtmp = _mm256_load_pd(_varX + j);
			DXtmp = _mm256_load_pd(_varDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);

			hDDX1 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX1 + j, hDDX1);
			tmp = _mm256_add_pd(_mm256_mul_pd(DXtmp, timeStep2), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, constantD2), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve2());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX2[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _hDDX1[j] * _timeStep4;
			//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

			Xtmp = _mm256_load_pd(_varX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX1 = _mm256_load_pd(_hDDX1 + j);

			hDDX2 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX2 + j, hDDX2);
			tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, timeStep4), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constantD2), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve3());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX3[j] = _varDDX[j] * _timeStep;
			//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			//_varDX[j] = _initDX[j] + _hDDX3[j];

			Xtmp = _mm256_load_pd(_initX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX2 = _mm256_load_pd(_hDDX2 + j);

			hDDX3 = _mm256_mul_pd(DDXtmp, timeStep);
			_mm256_store_pd(_hDDX3 + j, hDDX3);
			tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constantD2), DXtmp);
			tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			_mm256_store_pd(_varDX + j, _mm256_add_pd(DXtmp, hDDX3));
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve4());
		MeasuredRun(2, CalculateForces());

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//float sDDX = _hDDX2[j] + _hDDX3[j];
			//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

			Xtmp = _mm256_load_pd(_initX + j);
			DXtmp = _mm256_load_pd(_initDX + j);
			DDXtmp = _mm256_load_pd(_varDDX + j);
			hDDX1 = _mm256_load_pd(_hDDX1 + j);
			hDDX2 = _mm256_load_pd(_hDDX2 + j);
			hDDX3 = _mm256_load_pd(_hDDX3 + j);

			sDDX = _mm256_add_pd(hDDX2, hDDX3);
			tmp = _mm256_add_pd(_mm256_mul_pd(sDDX, constantD6), DXtmp);
			tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
			_mm256_store_pd(_varX + j, tmp);
			tmp = _mm256_add_pd(_mm256_add_pd(hDDX1, sDDX), sDDX);
			tmp = _mm256_add_pd(_mm256_mul_pd(DDXtmp, timeStep), tmp);
			tmp = _mm256_add_pd(_mm256_mul_pd(tmp, constantD6), DXtmp);
			_mm256_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve1());
		MeasuredRun(2, CalculateForces());

	}
	_testTimer.Stop(0);
#ifndef NOTIMER
	const int width = 16;
//	_testTimer.SetWidth(width);
	std::cout << "-----------------------------------\n";
	double t1 = _testTimer.Print(1, "Rotations: ");
	double t2 = _testTimer.Print(2, "Forces: ");
	double t3 = _testTimer.Print(3, "Integration: ");
	//_testTimer.Print(5, "Linksh:");
	//std::cout << std::setw(width) << "Summ: " << t1 + t2 + t3 << std::endl;
	_testTimer.Print(0, "Total: ");
#endif
}


void StressStrainCppIterativeSolverAVX::Solve(const int nIterations)
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	FTimer test_timer;
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

// virtual
void StressStrainCppIterativeSolverAVX::InitialSolve()
{
	MeasuredRun(1, _rotationSolver->InitialSolve());
	MeasuredRun(2, CalculateForces());
}


/**
* Расчет первой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolverAVX::Solve1()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	MeasuredRun(1, _rotationSolver->Solve1());
	MeasuredRun(2, CalculateForces());
}

/**
* Расчет второй стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolverAVX::Solve2()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep2;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep2 = _mm256_set1_pd(_timeStep2);
	constantD2 = _mm256_set1_pd(0.5);


	memcpy(_initX, _varX, sizeof(double)*_nVariables);
	memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX1[j] = _varDDX[j] * _timeStep;
		//_varX[j] += _varDX[j] * _timeStep2;
		//_varDX[j] += _hDDX1[j] * 0.5;

		__m256d Xtmp = _mm256_load_pd(_varX + j);
		__m256d DXtmp = _mm256_load_pd(_varDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);

		__m256d hDDX1 = _mm256_mul_pd(DDXtmp, timeStep);
		_mm256_store_pd(_hDDX1 + j, hDDX1);
		__m256d tmp = _mm256_add_pd(_mm256_mul_pd(DXtmp, timeStep2), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, constantD2), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve2());
	MeasuredRun(2, CalculateForces());

}

/**
* Расчет третьей стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolverAVX::Solve3()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d timeStep4;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	timeStep4 = _mm256_set1_pd(_timeStep4);
	constantD2 = _mm256_set1_pd(0.5);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX2[j] = _varDDX[j] * _timeStep;
		//_varX[j] += _hDDX1[j] * _timeStep4;
		//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

		__m256d Xtmp = _mm256_load_pd(_varX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX1 = _mm256_load_pd(_hDDX1 + j);
		__m256d hDDX2 = _mm256_mul_pd(DDXtmp, timeStep);

		_mm256_store_pd(_hDDX2 + j, hDDX2);
		__m256d tmp = _mm256_add_pd(_mm256_mul_pd(hDDX1, timeStep4), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constantD2), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve3());
	MeasuredRun(2, CalculateForces());
}

/**
* Расчет четвертой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolverAVX::Solve4()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d constantD2;

	timeStep = _mm256_set1_pd(_timeStep);
	constantD2 = _mm256_set1_pd(0.5);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//_hDDX3[j] = _varDDX[j] * _timeStep;
		//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
		//_varDX[j] = _initDX[j] + _hDDX3[j];

		__m256d Xtmp = _mm256_load_pd(_initX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX2 = _mm256_load_pd(_hDDX2 + j);

		__m256d hDDX3 = _mm256_mul_pd(DDXtmp, timeStep);
		_mm256_store_pd(_hDDX3 + j, hDDX3);
		__m256d tmp = _mm256_add_pd(_mm256_mul_pd(hDDX2, constantD2), DXtmp);
		tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		_mm256_store_pd(_varDX + j, _mm256_add_pd(DXtmp, hDDX3));
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve4());
	MeasuredRun(2, CalculateForces());
}

/**
* Расчет пятой стадии метода Рунге-Кутты
*/
// virtual
void StressStrainCppIterativeSolverAVX::Solve5()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

	__m256d timeStep;
	__m256d constantD6;

	timeStep = _mm256_set1_pd(_timeStep);
	constantD6 = _mm256_set1_pd(1 / 6.0);

	_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
	for (int j = 0; j < _nVariables; j += regSize)
	{
		//float sDDX = _hDDX2[j] + _hDDX3[j];
		//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
		//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

		__m256d Xtmp = _mm256_load_pd(_initX + j);
		__m256d DXtmp = _mm256_load_pd(_initDX + j);
		__m256d DDXtmp = _mm256_load_pd(_varDDX + j);
		__m256d hDDX1 = _mm256_load_pd(_hDDX1 + j);
		__m256d hDDX2 = _mm256_load_pd(_hDDX2 + j);
		__m256d hDDX3 = _mm256_load_pd(_hDDX3 + j);

		__m256d sDDX = _mm256_add_pd(hDDX2, hDDX3);
		__m256d tmp = _mm256_add_pd(_mm256_mul_pd(sDDX, constantD6), DXtmp);
		tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
		_mm256_store_pd(_varX + j, tmp);
		tmp = _mm256_add_pd(_mm256_add_pd(hDDX1, sDDX), sDDX);
		tmp = _mm256_add_pd(_mm256_mul_pd(DDXtmp, timeStep), tmp);
		tmp = _mm256_add_pd(_mm256_mul_pd(tmp, constantD6), DXtmp);
		_mm256_store_pd(_varDX + j, tmp);
	}
	_testTimer.Stop(3);

	MeasuredRun(1, _rotationSolver->Solve1());
	MeasuredRun(2, CalculateForces());
}

void StressStrainCppIterativeSolverAVX::CalculateForces()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	__declspec(align(32)) double strains[8], velocityStrains[8];

#pragma omp parallel for private (strains,velocityStrains) num_threads(_numThreads)
	for (int elementId1 = 0; elementId1 < _nElements; elementId1++)
	{
		
		double* accelerationVector = GetElementAcceleration(elementId1);

		memset(GetElementAcceleration(elementId1), 0u, sizeof(double)*vecStride2);
		memset(GetElementStress(elementId1), 0u, sizeof(double)*vecStride2);

		_testTimer.Start(5);

		// обход 6 связанных элементов x-,y-,z-,x+,y+,z+
		for (size_t side = 0; side < 6; side++)
		{
			size_t elementId2 = _linkedElements[6 * elementId1 + side];
			if (elementId2)
			{
				elementId2--;
				CalculateStrainsAVX(side, strains, velocityStrains, elementId1, elementId2);

				Vec3Ref linear_strains = MakeVec3(&strains[0]);
				Vec3Ref angular_strains = MakeVec3(&strains[0] + vecStride);
				Vec3Ref linear_vstrains = MakeVec3(&velocityStrains[0]);
				Vec3Ref angular_vstrains = MakeVec3(&velocityStrains[0] + vecStride);

				for (size_t component = 0; component < 3; component++)
				{
					GetElementStress(elementId1)[component] += strains[component];
					GetElementStress(elementId1)[component + vecStride] += strains[component + vecStride];
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
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
}
}

#endif