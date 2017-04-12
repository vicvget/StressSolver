#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>
#include <immintrin.h>

using std::string;
using std::vector;

namespace Stress
{

	class StressStrainCppIterativeSolverKNC:
		public StressStrainCppIterativeSolver
	{
		const int regSize = 8;

	public:

		// создает объект с заданными параметрами
		StressStrainCppIterativeSolverKNC
			(
			double* params,
			int* links,
			int nLinks,
			double *coordinates,
			int nElements,
			double gridStep,
			double timeStep,
			int numThreads,
			int stride
			);

		virtual
			~StressStrainCppIterativeSolverKNC();

#pragma region overriden

		virtual	void InitialSolve();

		virtual void Solve(const int nIteratons);
#ifndef DIRECT_INT
		virtual	void SolveFull(const int nIteratons);
#endif

		/**
		* Расчет первой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve1();

		/**
		* Расчет второй стадии метода Рунге-Кутты
		*/
		virtual
			void Solve2();

		/**
		* Расчет третьей стадии метода Рунге-Кутты
		*/
		virtual
			void Solve3();

		/**
		* Расчет четвертой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve4();

		/**
		* Расчет пятой стадии метода Рунге-Кутты
		*/
		virtual
			void Solve5();

#ifndef DIRECT_RHS
		virtual
			void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
#endif

#pragma endregion

		double df[12]; // debug
	};
#ifdef USE_KNC
	inline __m512d _mm512_loadu_pd(const double* a)
	{
		__m512d v_temp = _mm512_setzero_pd();
		v_temp = _mm512_loadunpacklo_pd(v_temp, &a[0]);
		v_temp = _mm512_loadunpackhi_pd(v_temp, &a[8]);

		return v_temp;
	}
#endif
}
