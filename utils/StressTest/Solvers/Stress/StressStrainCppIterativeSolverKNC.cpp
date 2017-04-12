#include "StressStrainCppIterativeSolverKNC.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;

namespace Stress
{

	StressStrainCppIterativeSolverKNC::StressStrainCppIterativeSolverKNC
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
		std::cout << "KNC SOLVER" << std::endl << std::flush;
	}

	// virtual
	StressStrainCppIterativeSolverKNC::~StressStrainCppIterativeSolverKNC()
	{
	}

	void StressStrainCppIterativeSolverKNC::Solve(const int nIterations)
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
			//std::cout << "Solve2 routine completed" << std::endl;
			Solve3();
			//std::cout << "Solve3 routine completed" << std::endl;
			Solve4();
			//std::cout << "Solve4 routine completed" << std::endl;
			Solve5();
			//std::cout << "Solve5 routine completed" << std::endl;
		}
	}


#ifndef DIRECT_INT
	// KNC-версия
	void StressStrainCppIterativeSolverKNC::SolveFull(const int nIterations)
	{
#if defined(USE_KNC) || defined(USE_KNL)

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		_iterationNumber = 0;
		_testTimer.Start(0);

		__m512d Xtmp, DXtmp, DDXtmp, tmp;
		__m512d hDDX1, hDDX2, hDDX3, sDDX, hpsDDX;

#ifndef OMP_SOLVE

		__m512d timeStep = _mm512_set1_pd(_timeStep);
		__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
		__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
		__m512d constantD2 = _mm512_set1_pd(0.5);
		__m512d constantD6 = _mm512_set1_pd(1 / 6.0);
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
				__m512d timeStep = _mm512_set1_pd(_timeStep);
				__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
				__m512d constantD2 = _mm512_set1_pd(0.5);
#endif
				Xtmp = _mm512_load_pd(_varX + j);
				DXtmp = _mm512_load_pd(_varDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);

				hDDX1 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX1 + j, hDDX1);
				tmp = _mm512_fmadd_pd(DXtmp, timeStep2, Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				tmp = _mm512_fmadd_pd(hDDX1, constantD2, DXtmp);
				_mm512_store_pd(_varDX + j, tmp);
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
				__m512d timeStep = _mm512_set1_pd(_timeStep);
				__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
				__m512d constantD2 = _mm512_set1_pd(0.5);
#endif
				Xtmp = _mm512_load_pd(_varX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX1 = _mm512_load_pd(_hDDX1 + j);

				hDDX2 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX2 + j, hDDX2);
				tmp = _mm512_fmadd_pd(hDDX1, timeStep4, Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				tmp = _mm512_fmadd_pd(hDDX2, constantD2, DXtmp);
				_mm512_store_pd(_varDX + j, tmp);
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
				__m512d timeStep = _mm512_set1_pd(_timeStep);
				__m512d constantD2 = _mm512_set1_pd(0.5);
#endif
				Xtmp = _mm512_load_pd(_initX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX2 = _mm512_load_pd(_hDDX2 + j);

				hDDX3 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX3 + j, hDDX3);
				tmp = _mm512_fmadd_pd(hDDX2, constantD2, DXtmp);
				tmp = _mm512_fmadd_pd(tmp, timeStep, Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				_mm512_store_pd(_varDX + j, _mm512_add_pd(DXtmp, hDDX3));
			}
			_testTimer.Stop(3);

			MeasuredRun(1, _rotationSolver->Solve4());
			MeasuredRun(2, CalculateForces());

			_testTimer.Start(3);
#ifdef OMP_SOLVE
#pragma omp parallel for num_threads(_numThreads) private(Xtmp, DXtmp, DDXtmp, hDDX1, hDDX2, hDDX3, sDDX, hpsDDX, tmp)
#endif
			for (int j = 0; j < _nVariables; j += regSize)
			{
				//float sDDX = _hDDX2[j] + _hDDX3[j];
				//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
				//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;
#ifdef OMP_SOLVE
				__m512d timeStep = _mm512_set1_pd(_timeStep);
				__m512d constantD6 = _mm512_set1_pd(1 / 6.0);
#endif
				Xtmp = _mm512_load_pd(_initX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX1 = _mm512_load_pd(_hDDX1 + j);
				hDDX2 = _mm512_load_pd(_hDDX2 + j);
				hDDX3 = _mm512_load_pd(_hDDX3 + j);

				sDDX = _mm512_add_pd(hDDX2, hDDX3);
				hpsDDX = _mm512_add_pd(hDDX1, sDDX);
				tmp = _mm512_fmadd_pd(hpsDDX, constantD6, DXtmp);
				tmp = _mm512_fmadd_pd(tmp, timeStep, Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				tmp = _mm512_add_pd(_mm512_add_pd(hDDX1, sDDX), sDDX);
				tmp = _mm512_fmadd_pd(DDXtmp, timeStep, tmp);
				tmp = _mm512_fmadd_pd(tmp, constantD6, DXtmp);
				_mm512_store_pd(_varDX + j, tmp);
			}
			_testTimer.Stop(3);

			MeasuredRun(1, _rotationSolver->Solve1());
			MeasuredRun(2, CalculateForces());
			CheckVelocitySumm();

		}
		_testTimer.Stop(0);

#endif
	}
#endif


	// virtual
	void StressStrainCppIterativeSolverKNC::InitialSolve()
	{
		MeasuredRun(1, _rotationSolver->InitialSolve());
		MeasuredRun(2, CalculateForces());
	}


	/**
	* Расчет первой стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve1()
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		MeasuredRun(1, _rotationSolver->Solve1());
		MeasuredRun(2, CalculateForces());
	}

	/**
	* Расчет второй стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve2()
	{
#if defined(USE_KNC) || defined(USE_KNL)

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

			__m512d timeStep = _mm512_set1_pd(_timeStep);
			__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
			__m512d constantD2 = _mm512_set1_pd(0.5);

			__m512d Xtmp = _mm512_load_pd(_varX + j);
			__m512d DXtmp = _mm512_load_pd(_varDX + j);
			__m512d DDXtmp = _mm512_load_pd(_varDDX + j);

			__m512d hDDX1 = _mm512_mul_pd(DDXtmp, timeStep);
			_mm512_store_pd(_hDDX1 + j, hDDX1);
			__m512d tmp = _mm512_add_pd(_mm512_mul_pd(DXtmp, timeStep2), Xtmp);
			_mm512_store_pd(_varX + j, tmp);
			tmp = _mm512_add_pd(_mm512_mul_pd(hDDX1, constantD2), DXtmp);
			_mm512_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve2());
		MeasuredRun(2, CalculateForces());
#endif
	}

	/**
	* Расчет третьей стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve3()
	{
#if defined(USE_KNC) || defined(USE_KNL)

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX2[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _hDDX1[j] * _timeStep4;
			//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

			__m512d timeStep = _mm512_set1_pd(_timeStep);
			__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
			__m512d constantD2 = _mm512_set1_pd(0.5);

			__m512d Xtmp = _mm512_load_pd(_varX + j);
			__m512d DXtmp = _mm512_load_pd(_initDX + j);
			__m512d DDXtmp = _mm512_load_pd(_varDDX + j);
			__m512d hDDX1 = _mm512_load_pd(_hDDX1 + j);
			__m512d hDDX2 = _mm512_mul_pd(DDXtmp, timeStep);

			_mm512_store_pd(_hDDX2 + j, hDDX2);
			__m512d tmp = _mm512_add_pd(_mm512_mul_pd(hDDX1, timeStep4), Xtmp);
			_mm512_store_pd(_varX + j, tmp);
			tmp = _mm512_add_pd(_mm512_mul_pd(hDDX2, constantD2), DXtmp);
			_mm512_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve3());
		MeasuredRun(2, CalculateForces());
#endif
	}

	/**
	* Расчет четвертой стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve4()
	{
#if defined(USE_KNC) || defined(USE_KNL)

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX3[j] = _varDDX[j] * _timeStep;
			//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			//_varDX[j] = _initDX[j] + _hDDX3[j];

			__m512d timeStep = _mm512_set1_pd(_timeStep);
			__m512d constantD2 = _mm512_set1_pd(0.5);

			__m512d Xtmp = _mm512_load_pd(_initX + j);
			__m512d DXtmp = _mm512_load_pd(_initDX + j);
			__m512d DDXtmp = _mm512_load_pd(_varDDX + j);
			__m512d hDDX2 = _mm512_load_pd(_hDDX2 + j);
			__m512d hDDX3 = _mm512_mul_pd(DDXtmp, timeStep);
			_mm512_store_pd(_hDDX3 + j, hDDX3);
			__m512d tmp = _mm512_add_pd(_mm512_mul_pd(hDDX2, constantD2), DXtmp);
			tmp = _mm512_add_pd(_mm512_mul_pd(tmp, timeStep), Xtmp);
			_mm512_store_pd(_varX + j, tmp);
			_mm512_store_pd(_varDX + j, _mm512_add_pd(DXtmp, hDDX3));
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve4());
		MeasuredRun(2, CalculateForces());
#endif
	}

	/**
	* Расчет пятой стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve5()
	{
#if defined(USE_KNC) || defined(USE_KNL)

		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//float sDDX = _hDDX2[j] + _hDDX3[j];
			//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

			__m512d timeStep = _mm512_set1_pd(_timeStep);
			__m512d constantD6 = _mm512_set1_pd(1 / 6.0);

			__m512d Xtmp = _mm512_load_pd(_initX + j);
			__m512d DXtmp = _mm512_load_pd(_initDX + j);
			__m512d DDXtmp = _mm512_load_pd(_varDDX + j);
			__m512d hDDX1 = _mm512_load_pd(_hDDX1 + j);
			__m512d hDDX2 = _mm512_load_pd(_hDDX2 + j);
			__m512d hDDX3 = _mm512_load_pd(_hDDX3 + j);
			__m512d sDDX = _mm512_add_pd(hDDX2, hDDX3);
			__m512d hpsDDX = _mm512_add_pd(hDDX1, sDDX);
			__m512d tmp = _mm512_add_pd(_mm512_mul_pd(hpsDDX, constantD6), DXtmp);
			tmp = _mm512_add_pd(_mm512_mul_pd(tmp, timeStep), Xtmp);
			_mm512_store_pd(_varX + j, tmp);
			tmp = _mm512_add_pd(_mm512_add_pd(hDDX1, sDDX), sDDX);
			tmp = _mm512_add_pd(_mm512_mul_pd(DDXtmp, timeStep), tmp);
			tmp = _mm512_add_pd(_mm512_mul_pd(tmp, constantD6), DXtmp);
			_mm512_store_pd(_varDX + j, tmp);
		}
		_testTimer.Stop(3);

		MeasuredRun(1, _rotationSolver->Solve1());
		MeasuredRun(2, CalculateForces());
#endif
	}

#ifndef DIRECT_RHS
	void StressStrainCppIterativeSolverKNC::CalculateStrains
	(
		size_t side,			// 0 = -x, 1 = x, 2 = -y, 3 = y, 4 = -z, 5 = z
		double *shiftStrains,		// выход деформаций
		double *velocityStrains,	// выход изм. скоростей
		size_t nodeId1,				// номер узла 1
		size_t nodeId2					// номер узла 2
	) const
	{
#if defined(USE_KNC) || defined(USE_KNL)

		// Start AVX code
		double* pmatA01 = GetRotationMatrix(nodeId1);
		double* pmatA02 = GetRotationMatrix(nodeId2);

		__m512d matA01row1;
		__m512d matA01row2;
		__m512d matA01row3;
		// vecStride must be 4
		if (nodeId1 % 2)
		{
			matA01row1 = _mm512_loadu_pd(pmatA01);
			matA01row2 = _mm512_load_pd(pmatA01 + vecStride);
			matA01row3 = _mm512_loadu_pd(pmatA01 + vecStride2);
		}
		else
		{
			matA01row1 = _mm512_load_pd(pmatA01);
			matA01row2 = _mm512_loadu_pd(pmatA01 + vecStride);
			matA01row3 = _mm512_load_pd(pmatA01 + vecStride2);
		}
		//std::cout << "Loaded matrix rows" << std::endl << std::flush;
		__m512d matA02el1;
		__m512d matA02el2;
		__m512d matA02el3;

		__m512d matA21row1;
		__m512d matA21row2;
		__m512d matA21row3;

		// matA02 column 1/3
		matA02el1 = _mm512_set1_pd(pmatA02[0]);
		matA02el2 = _mm512_set1_pd(pmatA02[vecStride]);
		matA02el3 = _mm512_set1_pd(pmatA02[vecStride2]);
		//std::cout << "Loaded matrix col elements" << std::endl << std::flush;

		// matA01.TMul(matA02)
		matA21row1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA21row1 = _mm512_fmadd_pd(matA01row2, matA02el2, matA21row1);
		matA21row1 = _mm512_fmadd_pd(matA01row3, matA02el3, matA21row1);

		//std::cout << "Multiply add" << std::endl << std::flush;

		// matA02 column 2/3
		matA02el1 = _mm512_set1_pd(pmatA02[1]);
		matA02el2 = _mm512_set1_pd(pmatA02[vecStride + 1]);
		matA02el3 = _mm512_set1_pd(pmatA02[vecStride2 + 1]);

		// matA01.TMul(matA02)
		matA21row2 = _mm512_mul_pd(matA01row1, matA02el1);
		matA21row2 = _mm512_fmadd_pd(matA01row2, matA02el2, matA21row2);
		matA21row2 = _mm512_fmadd_pd(matA01row3, matA02el3, matA21row2);

		// matA02 column 3/3
		matA02el1 = _mm512_set1_pd(pmatA02[2]);
		matA02el2 = _mm512_set1_pd(pmatA02[vecStride + 2]);
		matA02el3 = _mm512_set1_pd(pmatA02[vecStride2 + 2]);

		// matA01.TMul(matA02)
		matA21row3 = _mm512_mul_pd(matA01row1, matA02el1);
		matA21row3 = _mm512_fmadd_pd(matA01row2, matA02el2, matA21row3);
		matA21row3 = _mm512_fmadd_pd(matA01row3, matA02el3, matA21row3);

		// Матрица A_{21} сформирована

		double* vecC1 = GetRadiusVector(side);
		__m512d ivecC1 = _mm512_loadu_pd(vecC1);
		__m512d vecDP = _mm512_sub_pd(
			_mm512_load_pd(GetElementShift(nodeId1)),
			_mm512_load_pd(GetElementShift(nodeId2))); // P1-P2

		matA02el1 = _mm512_set1_pd(vecC1[0]);
		matA02el2 = _mm512_set1_pd(vecC1[1]);
		matA02el3 = _mm512_set1_pd(vecC1[2]);

		__m512d res;
		res = _mm512_fmadd_pd(matA21row1, matA02el1, ivecC1);
		res = _mm512_fmadd_pd(matA21row2, matA02el2, res);
		res = _mm512_fmadd_pd(matA21row3, matA02el3, res);

		__declspec(align(64)) double tmp[8] = { 0 };
		_mm512_store_pd(tmp, vecDP);
		matA02el1 = _mm512_set1_pd(tmp[0]);
		matA02el2 = _mm512_set1_pd(tmp[1]);
		matA02el3 = _mm512_set1_pd(tmp[2]);

		res = _mm512_fmadd_pd(matA01row1, matA02el1, res);
		res = _mm512_fmadd_pd(matA01row2, matA02el2, res);
		res = _mm512_fmadd_pd(matA01row3, matA02el3, res);

		_mm512_store_pd(shiftStrains, res); // получено SL, линейные компоненты

		// Расчет VL
		// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
		// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
		// vecC2 = -vecC1
		// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

		__declspec(align(64)) double cp1[8] = { 0 };
		__declspec(align(64)) double cp2[8] = { 0 };

		CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
		CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]

		res = _mm512_load_pd(&cp1[0]);
		__m512d vecDV = _mm512_sub_pd(
			_mm512_load_pd(GetElementVelocity(nodeId1)),
			_mm512_load_pd(GetElementVelocity(nodeId2))); // V1-V2

		matA02el1 = _mm512_set1_pd(cp2[0]);
		matA02el2 = _mm512_set1_pd(cp2[1]);
		matA02el3 = _mm512_set1_pd(cp2[2]);

		res = _mm512_fmadd_pd(matA21row1, matA02el1, res);
		res = _mm512_fmadd_pd(matA21row2, matA02el2, res);
		res = _mm512_fmadd_pd(matA21row3, matA02el3, res);

		_mm512_store_pd(tmp, vecDV);
		matA02el1 = _mm512_set1_pd(tmp[0]);
		matA02el2 = _mm512_set1_pd(tmp[1]);
		matA02el3 = _mm512_set1_pd(tmp[2]);

		res = _mm512_fmadd_pd(matA01row1, matA02el1, res);
		res = _mm512_fmadd_pd(matA01row2, matA02el2, res);
		res = _mm512_fmadd_pd(matA01row3, matA02el3, res);

		_mm512_store_pd(velocityStrains, res); // получено VL, линейные компоненты

		double* sp1 = GetElementShiftAngular(nodeId1);
		double* sp2 = GetElementShiftAngular(nodeId2);
		double* vp1 = GetElementVelocityAngular(nodeId1);
		double* vp2 = GetElementVelocityAngular(nodeId2);
		for (size_t i = 0; i < 3; i++)
		{
			shiftStrains[i + vecStride] = sp1[i] - sp2[i];
			velocityStrains[i + vecStride] = vp1[i] - vp2[i];
		}

#endif
	}
#endif
}
