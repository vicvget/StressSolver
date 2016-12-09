#ifdef USE_KNC

#include "StressStrainCppIterativeSolverKNC.h"
#include "../../AdditionalModules/fmath/Matrix3x3.h"
#include "../../AdditionalModules/fmath/Matrix3x4.h"

using MathHelpers::Vec3;
using MathHelpers::Vec3Ref;
using MathHelpers::MakeVec3;
using MathHelpers::Mat3;
using MathHelpers::Mat3x4;


//__m512i izmm00 = _mm512_load_epi64(pmtx);    //A2A1
//__m512i izmm01 = _mm512_load_epi64(pmtx+8);	 //B1A3
//__m512i izmm02 = _mm512_load_epi64(pmtx+16); //B3B2
//__m512i izmm03 = _mm512_alignr_epi32(izmm00, izmm01, 8); //A1B1
//__m512i izmm04 = _mm512_alignr_epi32(izmm01, izmm02, 8); //A3B3
//
//__m512d zmm00 = (__m512d)_mm512_alignr_epi32(izmm03, izmm03, 8); //B1A1
//__m512d zmm01 = (__m512d)_mm512_alignr_epi32(izmm02, izmm00, 8); //B2A2
//__m512d zmm02 = (__m512d)_mm512_alignr_epi32(izmm04, izmm04, 8); //B3A3
//
////Вместо D1C1,D2C2,D3C3
////нам надо получить С1B1,C2B2,C3B3
//
//izmm03 = _mm512_load_epi64(pmtr);    //C2C1
//izmm04 = _mm512_load_epi64(pmtr+8);  //D1C3 // D1 не используется
//izmm05 = _mm512_alignr_epi32(izmm02, izmm03, 8); //B2C2
//
//zmm03 = (__m512d)_mm512_alignr_epi32(izmm03, izmm01, 8); //C1B1
//zmm04 = (__m512d)_mm512_alignr_epi32(izmm02, izmm02, 8); //C2B2
//zmm05 = (__m512d)_mm512_alignr_epi32(izmm04, izmm02, 8); //C3B3
//
//// остальная часть кода такая же

namespace Stress
{

	inline __m512d _mm512_loadu_pd(const double* a)
	{
		__m512d v_temp = _mm512_setzero_pd();
		v_temp = _mm512_loadunpacklo_pd(v_temp, &a[0]);
		v_temp = _mm512_loadunpackhi_pd(v_temp, &a[8]);

		return v_temp;
	}

	void StressStrainCppIterativeSolverKNC::CalculateStrains
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
		matA02el1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA01row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA01row3, matA02el3);

		matA21row1 = _mm512_add_pd(matA02el1, _mm512_add_pd(matA02el2, matA02el3));

		//std::cout << "Multiply add" << std::endl << std::flush;

		// matA02 column 2/3
		matA02el1 = _mm512_set1_pd(pmatA02[1]);
		matA02el2 = _mm512_set1_pd(pmatA02[vecStride + 1]);
		matA02el3 = _mm512_set1_pd(pmatA02[vecStride2 + 1]);

		// matA01.TMul(matA02)
		matA02el1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA01row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA01row3, matA02el3);

		matA21row2 = _mm512_add_pd(matA02el1, _mm512_add_pd(matA02el2, matA02el3));

		// matA02 column 3/3
		matA02el1 = _mm512_set1_pd(pmatA02[2]);
		matA02el2 = _mm512_set1_pd(pmatA02[vecStride + 2]);
		matA02el3 = _mm512_set1_pd(pmatA02[vecStride2 + 2]);

		// matA01.TMul(matA02)
		matA02el1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA01row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA01row3, matA02el3);

		matA21row3 = _mm512_add_pd(matA02el1, _mm512_add_pd(matA02el2, matA02el3));
		// Матрица A_{21} сформирована

		//std::cout << "Multiply add completed" << std::endl << std::flush;

		double* vecC1 = GetRadiusVector(side);
		__m512d ivecC1 = _mm512_loadu_pd(vecC1);
		__m512d vecDP = _mm512_sub_pd(
			_mm512_load_pd(GetElementShift(nodeId1)),
			_mm512_load_pd(GetElementShift(nodeId2))); // P1-P2
		//std::cout << "Radius vector loaded" << std::endl << std::flush;


		//__declspec(align(64)) double tmp[8]
		//double tmp[8]  __attribute__((aligned(64)));
		double* tmp = _buffer;
		_mm512_store_pd(tmp, ivecC1);
		//std::cout << "Stored to tmp" << std::endl << std::flush;

		matA02el1 = _mm512_set1_pd(vecC1[0]);
		matA02el2 = _mm512_set1_pd(vecC1[1]);
		matA02el3 = _mm512_set1_pd(vecC1[2]);

		matA02el1 = _mm512_mul_pd(matA21row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA21row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA21row3, matA02el3);

		__m512d mul1 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

		_mm512_store_pd(tmp, vecDP);
		matA02el1 = _mm512_set1_pd(tmp[0]);
		matA02el2 = _mm512_set1_pd(tmp[1]);
		matA02el3 = _mm512_set1_pd(tmp[2]);

		matA02el1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA01row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA01row3, matA02el3);

		__m512d mul2 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);
		__m512d res = _mm512_add_pd(_mm512_add_pd(ivecC1, mul1), mul2);

		_mm512_store_pd(shiftStrains, res); // получено SL, линейные компоненты

		// Расчет VL
		// переводим вектор разницы линейных скоростей точек связи С2-С1 в СК1
		// Vec3 VecT1 = matA01.Tmul(vecV1 - vecV2) + vecW1.Cross(vecC1) - matA12*(vecW2.Cross(vecC2));
		// vecC2 = -vecC1
		// vecC2.Cross(vecW1) = -vecC1.Cross(vecW1)

		//__declspec(align(64)) double cp1[8];
		//__declspec(align(64)) double cp2[8];

		//double cp1[8] __attribute__((aligned(64)));
		//double cp2[8] __attribute__((aligned(64)));
		double* cp1 = _buffer + 8;
		double* cp2 = _buffer + 16;
		CrossProduct(GetElementVelocityAngular(nodeId1), GetRadiusVector(side), cp1);	// [w1 x c1]
		CrossProduct(GetElementVelocityAngular(nodeId2), GetRadiusVector(side), cp2);	// -[w2 x c2] = [w2 x c1]

		__m512d cp1r = _mm512_load_pd(&cp1[0]);
		__m512d vecDV = _mm512_sub_pd(
			_mm512_load_pd(GetElementVelocity(nodeId1)),
			_mm512_load_pd(GetElementVelocity(nodeId2))); // V1-V2

		matA02el1 = _mm512_set1_pd(cp2[0]);
		matA02el2 = _mm512_set1_pd(cp2[1]);
		matA02el3 = _mm512_set1_pd(cp2[2]);

		matA02el1 = _mm512_mul_pd(matA21row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA21row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA21row3, matA02el3);

		mul1 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

		_mm512_store_pd(tmp, vecDV);
		matA02el1 = _mm512_set1_pd(tmp[0]);
		matA02el2 = _mm512_set1_pd(tmp[1]);
		matA02el3 = _mm512_set1_pd(tmp[2]);

		matA02el1 = _mm512_mul_pd(matA01row1, matA02el1);
		matA02el2 = _mm512_mul_pd(matA01row2, matA02el2);
		matA02el3 = _mm512_mul_pd(matA01row3, matA02el3);

		mul2 = _mm512_add_pd(_mm512_add_pd(matA02el1, matA02el2), matA02el3);

		res = _mm512_add_pd(_mm512_add_pd(cp1r, mul1), mul2);

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
	}

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

	// AVX-версия
	void StressStrainCppIterativeSolverKNC::SolveFull(const int nIterations)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		_iterationNumber = 0;
		_testTimer.Start(0);

		__m512d Xtmp, DXtmp, DDXtmp, tmp;
		__m512d hDDX1, hDDX2, hDDX3, sDDX;

		//__declspec(align(64)) double buffer[8];
		//buffer[0] = _timeStep;
		//buffer[1] = _timeStep2;
		//buffer[2] = _timeStep4;
		//buffer[3] = 0.5;
		//buffer[4] = 1 / 6.0;

		//__m512d timeStep = _mm512_extload_pd(&buffer[0],_MM_UPCONV_PD_NONE,_MM_BROADCAST_1X8,0);
		//__m512d timeStep2 = _mm512_extload_pd(&buffer[1],_MM_UPCONV_PD_NONE,_MM_BROADCAST_1X8,0);
		//__m512d timeStep4 = _mm512_extload_pd(&buffer[2],_MM_UPCONV_PD_NONE,_MM_BROADCAST_1X8,0);
		//__m512d constantD2 = _mm512_extload_pd(&buffer[3],_MM_UPCONV_PD_NONE,_MM_BROADCAST_1X8,0);
		//__m512d constantD6 = _mm512_extload_pd(&buffer[4],_MM_UPCONV_PD_NONE,_MM_BROADCAST_1X8,0);

		__m512d timeStep = _mm512_set1_pd(_buffer[8 * 3 + 0]);
		__m512d timeStep2 = _mm512_set1_pd(_buffer[8 * 3 + 1]);
		__m512d timeStep4 = _mm512_set1_pd(_buffer[8 * 3 + 2]);
		__m512d constantD2 = _mm512_set1_pd(_buffer[8 * 3 + 3]);
		__m512d constantD6 = _mm512_set1_pd(_buffer[8 * 3 + 4]);

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

				Xtmp = _mm512_load_pd(_varX + j);
				DXtmp = _mm512_load_pd(_varDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);

				hDDX1 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX1 + j, hDDX1);
				tmp = _mm512_add_pd(_mm512_mul_pd(DXtmp, timeStep2), Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				tmp = _mm512_add_pd(_mm512_mul_pd(hDDX1, constantD2), DXtmp);
				_mm512_store_pd(_varDX + j, tmp);
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

				Xtmp = _mm512_load_pd(_varX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX1 = _mm512_load_pd(_hDDX1 + j);

				hDDX2 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX2 + j, hDDX2);
				tmp = _mm512_add_pd(_mm512_mul_pd(hDDX1, timeStep4), Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				tmp = _mm512_add_pd(_mm512_mul_pd(hDDX2, constantD2), DXtmp);
				_mm512_store_pd(_varDX + j, tmp);
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

				Xtmp = _mm512_load_pd(_initX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX2 = _mm512_load_pd(_hDDX2 + j);

				hDDX3 = _mm512_mul_pd(DDXtmp, timeStep);
				_mm512_store_pd(_hDDX3 + j, hDDX3);
				tmp = _mm512_add_pd(_mm512_mul_pd(hDDX2, constantD2), DXtmp);
				tmp = _mm512_add_pd(_mm512_mul_pd(tmp, timeStep), Xtmp);
				_mm512_store_pd(_varX + j, tmp);
				_mm512_store_pd(_varDX + j, _mm512_add_pd(DXtmp, hDDX3));
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

				Xtmp = _mm512_load_pd(_initX + j);
				DXtmp = _mm512_load_pd(_initDX + j);
				DDXtmp = _mm512_load_pd(_varDDX + j);
				hDDX1 = _mm512_load_pd(_hDDX1 + j);
				hDDX2 = _mm512_load_pd(_hDDX2 + j);
				hDDX3 = _mm512_load_pd(_hDDX3 + j);

				sDDX = _mm512_add_pd(hDDX2, hDDX3);
				tmp = _mm512_add_pd(_mm512_mul_pd(sDDX, constantD6), DXtmp);
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


	void StressStrainCppIterativeSolverKNC::Solve(const int nIterations)
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		_iterationNumber = 0;
		_testTimer.Start(0);
		//std::cout << "Solve routine entered" << std::endl;
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

	// virtual
	void StressStrainCppIterativeSolverKNC::InitialSolve()
	{
		*(_buffer + 8 * 3 + 0) = _timeStep;
		*(_buffer + 8 * 3 + 1) = _timeStep2;
		*(_buffer + 8 * 3 + 2) = _timeStep4;
		*(_buffer + 8 * 3 + 3) = 0.5;
		*(_buffer + 8 * 3 + 4) = 1 / 6.0;
		//timeStep = _mm512_set1_pd(_buffer[8*3+0]);
		//timeStep2 = _mm512_set1_pd(_buffer[8*3+1]);
		//timeStep4 = _mm512_set1_pd(_buffer[8*3+2]);
		//constantD2 = _mm512_set1_pd(_buffer[8*3+3]);
		//constantD6 = _mm512_set1_pd(_buffer[8*3+4]);

		MeasuredRun(1, _rotationSolver->InitialSolve());
		std::cout << "ISolve completed" << std::endl << std::flush;
		MeasuredRun(2, CalculateForces());
		std::cout << "CForces completed" << std::endl << std::flush;
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
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		__m512d timeStep = _mm512_set1_pd(_timeStep);
		__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
		__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
		__m512d constantD2 = _mm512_set1_pd(0.5);
		__m512d constantD6 = _mm512_set1_pd(1 / 6.0);

		memcpy(_initX, _varX, sizeof(double)*_nVariables);
		memcpy(_initDX, _varDX, sizeof(double)*_nVariables);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX1[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _varDX[j] * _timeStep2;
			//_varDX[j] += _hDDX1[j] * 0.5;

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

	}

	/**
	* Расчет третьей стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve3()
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		__m512d timeStep = _mm512_set1_pd(_timeStep);
		__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
		__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
		__m512d constantD2 = _mm512_set1_pd(0.5);
		__m512d constantD6 = _mm512_set1_pd(1 / 6.0);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX2[j] = _varDDX[j] * _timeStep;
			//_varX[j] += _hDDX1[j] * _timeStep4;
			//_varDX[j] = _initDX[j] + _hDDX2[j] * 0.5;

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
	}

	/**
	* Расчет четвертой стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve4()
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		__m512d timeStep = _mm512_set1_pd(_timeStep);
		__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
		__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
		__m512d constantD2 = _mm512_set1_pd(0.5);
		__m512d constantD6 = _mm512_set1_pd(1 / 6.0);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//_hDDX3[j] = _varDDX[j] * _timeStep;
			//_varX[j] = _initX[j] + (_initDX[j] + _hDDX2[j] * 0.5) * _timeStep;
			//_varDX[j] = _initDX[j] + _hDDX3[j];

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
	}

	/**
	* Расчет пятой стадии метода Рунге-Кутты
	*/
	// virtual
	void StressStrainCppIterativeSolverKNC::Solve5()
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		__m512d timeStep = _mm512_set1_pd(_timeStep);
		__m512d timeStep2 = _mm512_set1_pd(_timeStep2);
		__m512d timeStep4 = _mm512_set1_pd(_timeStep4);
		__m512d constantD2 = _mm512_set1_pd(0.5);
		__m512d constantD6 = _mm512_set1_pd(1 / 6.0);

		_testTimer.Start(3);
#pragma omp parallel for num_threads(_numThreads)
		for (int j = 0; j < _nVariables; j += regSize)
		{
			//float sDDX = _hDDX2[j] + _hDDX3[j];
			//_varX[j] = _initX[j] + (_initDX[j] + sDDX / 6.0) * _timeStep;
			//_varDX[j] = _initDX[j] + (_hDDX1[j] + sDDX + sDDX + _varDDX[j] * _timeStep) / 6.0;

			__m512d Xtmp = _mm512_load_pd(_initX + j);
			__m512d DXtmp = _mm512_load_pd(_initDX + j);
			__m512d DDXtmp = _mm512_load_pd(_varDDX + j);
			__m512d hDDX1 = _mm512_load_pd(_hDDX1 + j);
			__m512d hDDX2 = _mm512_load_pd(_hDDX2 + j);
			__m512d hDDX3 = _mm512_load_pd(_hDDX3 + j);
			__m512d sDDX = _mm512_add_pd(hDDX2, hDDX3);
			__m512d tmp = _mm512_add_pd(_mm512_mul_pd(sDDX, constantD6), DXtmp);
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
	}

	void StressStrainCppIterativeSolverKNC::CalculateForces()
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		static int it = 0;

		//__declspec(align(64)) double strains[8];
		//__declspec(align(64)) double velocityStrains[8];
		double velocityStrains[8] __attribute__((aligned(64)));
		double strains[8] __attribute__((aligned(64)));
		for (size_t elementId1 = 0; elementId1 < _nElements; elementId1++)
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
					CalculateStrains(side, strains, velocityStrains, elementId1, elementId2);

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
					Vec3 force = -linear_vstrains * _dampingFactorLinear - linear_strains * _elasticFactorLinear;
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