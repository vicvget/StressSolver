#ifndef USE_KNC
//#define OMP_SOLVE
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

	void StressStrainCppIterativeSolverAVX::Solve(const int nIterations)
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

			Solve2();
			Solve3();
			Solve4();
			Solve5();
			CheckVelocitySumm();
		}
	}

#ifndef DIRECT_INT
	// AVX-версия
	void StressStrainCppIterativeSolverAVX::SolveFull(const int nIterations)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);


		_iterationNumber = 0;
		_testTimer.Start(0);

		__m256d Xtmp, DXtmp, DDXtmp, tmp;
		__m256d hDDX1, hDDX2, hDDX3, sDDX, hpsDDX1;

#ifndef OMP_SOLVE
		__m256d timeStep = _mm256_set1_pd(_timeStep);
		__m256d timeStep2 = _mm256_set1_pd(_timeStep2);
		__m256d timeStep4 = _mm256_set1_pd(_timeStep4);
		__m256d constantD2 = _mm256_set1_pd(0.5);
		__m256d constantD6 = _mm256_set1_pd(1 / 6.0);
#endif

		while (_iterationNumber != nIterations && _rotationSolver->IsValid())
		{
			_iterationNumber++;
			_nIteration++;

			_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
			memcpy(_initX, _varX, sizeof(double)*_nVariables);
			memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

			_testTimer.Start(3);
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)  private(Xtmp, DXtmp, DDXtmp, hDDX1, tmp)
#endif			
			for (int j = 0; j < _nVariables; j += regSize)
			{
				//_hDDX1[j] = _varDDX[j] * _timeStep;
				//_varX[j] += _varDX[j] * _timeStep2;
				//_varDX[j] += _hDDX1[j] * 0.5;

#ifdef OMP_SOLVE
				__m256d timeStep = _mm256_set1_pd(_timeStep);
				__m256d timeStep2 = _mm256_set1_pd(_timeStep2);
				__m256d constantD2 = _mm256_set1_pd(0.5);
#endif
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
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)  private(Xtmp, DXtmp, DDXtmp, hDDX1, hDDX2, tmp)
#endif
			for (int j = 0; j < _nVariables; j += regSize)
			{
				//_hDDX2[j] = _varDDX[j] * _timeStep;
				//_varX[j] += _hDDX1[j] * _timeStep4;
				//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

#ifdef OMP_SOLVE
				__m256d timeStep = _mm256_set1_pd(_timeStep);
				__m256d timeStep4 = _mm256_set1_pd(_timeStep4);
				__m256d constantD2 = _mm256_set1_pd(0.5);
#endif

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
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads) private(Xtmp, DXtmp, DDXtmp, hDDX3, hDDX2, tmp)
#endif
			for (int j = 0; j < _nVariables; j += regSize)
			{
				//_hDDX3[j] = _varDDX[j] * _timeStep;
				//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
				//_varDX[j] = _initDX[j] + _hDDX3[j];

#ifdef OMP_SOLVE
				__m256d timeStep = _mm256_set1_pd(_timeStep);
				__m256d constantD2 = _mm256_set1_pd(0.5);
#endif

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
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads) private(Xtmp, DXtmp, DDXtmp, hDDX1, hDDX2, hDDX3, sDDX, hpsDDX1, tmp)
#endif
			for (int j = 0; j < _nVariables; j += regSize)
			{
				//float sDDX = _hDDX2[j] + _hDDX3[j];
				//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
				//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

#ifdef OMP_SOLVE
				__m256d timeStep = _mm256_set1_pd(_timeStep);
				__m256d constantD6 = _mm256_set1_pd(1 / 6.0);
#endif

				Xtmp = _mm256_load_pd(_initX + j);
				DXtmp = _mm256_load_pd(_initDX + j);
				DDXtmp = _mm256_load_pd(_varDDX + j);
				hDDX1 = _mm256_load_pd(_hDDX1 + j);

				hDDX2 = _mm256_load_pd(_hDDX2 + j);
				hDDX3 = _mm256_load_pd(_hDDX3 + j);
				sDDX = _mm256_add_pd(hDDX2, hDDX3);
				hpsDDX1 = _mm256_add_pd(hDDX1, sDDX);

				tmp = _mm256_add_pd(_mm256_mul_pd(hpsDDX1, constantD6), DXtmp);
				tmp = _mm256_add_pd(_mm256_mul_pd(tmp, timeStep), Xtmp);
				_mm256_store_pd(_varX + j, tmp);

				tmp = _mm256_add_pd(hpsDDX1, sDDX);
				tmp = _mm256_add_pd(_mm256_mul_pd(DDXtmp, timeStep), tmp);
				tmp = _mm256_add_pd(_mm256_mul_pd(tmp, constantD6), DXtmp);
				_mm256_store_pd(_varDX + j, tmp);
			}
			_testTimer.Stop(3);

			MeasuredRun(1, _rotationSolver->Solve1());
			MeasuredRun(2, CalculateForces());

			CheckVelocitySumm();
		}
		_testTimer.Stop(0);
	}
#endif

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

		_testTimer.Start(3);


		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX1[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _varDX[j] * _timeStep2;
			//_varDX[j] += _hDDX1[j] * 0.5;

			__m256d timeStep = _mm256_set1_pd(_timeStep);
			__m256d timeStep2 = _mm256_set1_pd(_timeStep2);
			__m256d constantD2 = _mm256_set1_pd(0.5);

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

		_testTimer.Start(3);


#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX2[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _hDDX1[j] * _timeStep4;
			//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;
			__m256d timeStep = _mm256_set1_pd(_timeStep);
			__m256d timeStep4 = _mm256_set1_pd(_timeStep4);
			__m256d constantD2 = _mm256_set1_pd(0.5);

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

		_testTimer.Start(3);


#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX3[j] = _varDDX[j] * _timeStep;
			//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			//_varDX[j] = _initDX[j] + _hDDX3[j];

			__m256d timeStep = _mm256_set1_pd(_timeStep);
			__m256d constantD2 = _mm256_set1_pd(0.5);
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

		_testTimer.Start(3);


#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads)
#endif
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//float sDDX = _hDDX2[j] + _hDDX3[j];
			//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
			__m256d timeStep = _mm256_set1_pd(_timeStep);
			__m256d constantD6 = _mm256_set1_pd(1 / 6.0);

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

	void StressStrainCppIterativeSolverAVX::CalculateStrains
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

		double* rv = GetRadiusVector(side);
		__m256d ivecC1 = _mm256_load_pd(rv);
		__m256d vecDP = _mm256_sub_pd(
			_mm256_load_pd(GetElementShift(nodeId1)),
			_mm256_load_pd(GetElementShift(nodeId2))); // P1-P2

		matA02el1 = _mm256_set1_pd(rv[0]);
		matA02el2 = _mm256_set1_pd(rv[1]);
		matA02el3 = _mm256_set1_pd(rv[2]);

		matA02el1 = _mm256_mul_pd(matA21row1, matA02el1);
		matA02el2 = _mm256_mul_pd(matA21row2, matA02el2);
		matA02el3 = _mm256_mul_pd(matA21row3, matA02el3);

		__m256d mul1 = _mm256_add_pd(_mm256_add_pd(matA02el1, matA02el2), matA02el3);
		__declspec(align(32)) double tmp[4];
		_mm256_store_pd(tmp, vecDP);
		matA02el1 = _mm256_set1_pd(tmp[0]);
		matA02el2 = _mm256_set1_pd(tmp[1]);
		matA02el3 = _mm256_set1_pd(tmp[2]);

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

		_mm256_store_pd(tmp, vecDV);
		matA02el1 = _mm256_set1_pd(tmp[0]);
		matA02el2 = _mm256_set1_pd(tmp[1]);
		matA02el3 = _mm256_set1_pd(tmp[2]);

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
    
    void StressStrainCppIterativeSolverAVX::CalculateForces()
    {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	static int it = 0;

	__declspec(align(64)) double strains[8], velocityStrains[8],force0[4],force1[4], force2[4],torque[4];


	memset(_dataInternal + (_nElements * vecStride2 * 2), 0u, sizeof(double)*vecStride2*_nElements);
	memset(_stress, 0u, sizeof(double)*vecStride2*_nElements);
//#define timeMeas

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
#ifdef timeMeas
                size_t timer1 = __rdtsc();
#endif

				CalculateStrains(dof, strains, velocityStrains, elementId1, elementId2);
#ifdef timeMeas
                size_t timer2 = __rdtsc();
#endif



                double * stressFactors = _elementStressFactorCache + (elementId1 * vecStride);
                double * elementStress = _stress + (elementId1 * vecStride2);
                double * elementStressAngularId1 = _stress + (elementId1 * vecStride2 + vecStride);
                double * elementStressAngularId2 = _stress + (elementId2 * vecStride2 + vecStride);

				// нормальные напряжения
				elementStress[dof] += strains[dof] * stressFactors[dof] * _stressScalingFactors[dof];
				elementStress[dof] += strains[dof] * stressFactors[dof] * _stressScalingFactors[dof];

				// степени свободы смещений, участвующих в создании касательных напряжений
				int dof0 = exclusive_dofs[dof][0];
				int dof1 = exclusive_dofs[dof][1];

				elementStressAngularId1[dof0] += strains[dof1] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];
				elementStressAngularId1[dof1] += strains[dof0] * GetElementStressFactors(elementId1)[dof] * _stressScalingFactors[dof];

				elementStressAngularId2[dof0] += strains[dof1] * (_elementStressFactorCache + (elementId2 * vecStride))[dof] * _stressScalingFactors[dof];
				elementStressAngularId2[dof1] += strains[dof0] * (_elementStressFactorCache + (elementId2 * vecStride))[dof] * _stressScalingFactors[dof];

				// сила и момент из полученных деформаций
                double dfl = -_dampingFactorLinear;
                double dfa = -_dampingFactorAngular;

                vecDoubleMulSub(velocityStrains, strains, _dampingFactorLinear, _elasticFactorLinear, force1);
                vecDoubleMulSub(velocityStrains + vecStride, strains + vecStride, _dampingFactorAngular, _elasticFactorAngular, torque);
                
                


                
                double * matA01 = _dataRotationMtx + (elementId1 * matStride);
                double * matA02 = _dataRotationMtx + (elementId2 * matStride);
                double * accel1 = _dataInternal + (_nElements * vecStride2 * 2 + elementId1 * vecStride2);
                double * accel2 = _dataInternal + (_nElements * vecStride2 * 2 + elementId2 * vecStride2);
                matVecMul3x4(matA01,force1,force0);
                matVecTMul3x4(matA02, force0, force2);
                
                for (size_t i = 0; i < 3; i++)
                {
                    accel1[i] += force0[i];
                    accel2[i] -= force0[i];
                }

                
                
                double * vR = _radiusVectors + dof * vecStride;
                double * AccelAngul1 = _dataInternal + (_nElements * vecStride2 * 2 + elementId1 * vecStride2 + vecStride);
                double * AccelAngul2 = _dataInternal + (_nElements * vecStride2 * 2 + elementId2 * vecStride2 + vecStride);
                
                AccelAngul1[0] +=  vR[1] * force1[2] - vR[2] * force1[1] + torque[0];
                AccelAngul1[1] += -vR[0] * force1[2] + vR[2] * force1[0] + torque[1];
                AccelAngul1[2] +=  vR[0] * force1[1] - vR[1] * force1[0] + torque[2];

                AccelAngul2[0] +=  vR[1] * force2[2] - vR[2] * force2[1] - torque[0];
                AccelAngul2[1] += -vR[0] * force2[2] + vR[2] * force2[0] - torque[1];
                AccelAngul2[2] +=  vR[0] * force2[1] - vR[1] * force2[0] - torque[2];
                
#ifdef timeMeas
                size_t timer3 = __rdtsc();
                std::cout << "Strains time\t" << timer2 - timer1 << std::endl;
                std::cout << "Force time\t"   << timer3 - timer2 << std::endl;
#endif
			}
		}
	}
	ApplyBoundary(); // модифицирует силы и моменты
	ApplyMass();	 // вычисляет ускорения делением сил на массы и моментов на моменты инерции
    }
}
//#undef OMP_SOLVE

#endif