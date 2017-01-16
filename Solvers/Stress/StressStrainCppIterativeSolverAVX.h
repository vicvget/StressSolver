#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>
#include <immintrin.h>

using std::string;
using std::vector;

namespace Stress
{

	class StressStrainCppIterativeSolverAVX:
		public StressStrainCppIterativeSolver
	{
		const int regSize = 4;

	public:

		// ������� ������ � ��������� �����������
		StressStrainCppIterativeSolverAVX
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
			~StressStrainCppIterativeSolverAVX();

#pragma region overriden

		virtual	void InitialSolve();

		virtual void Solve(const int nIteratons);
#ifndef DIRECT_INT
		//virtual	void SolveFull(const int nIteratons);
#endif
		virtual	void SolveFull(const int nIteratons);
		/**
		* ������ ������ ������ ������ �����-�����
		*/
		virtual
			void Solve1();

		/**
		* ������ ������ ������ ������ �����-�����
		*/
		virtual
			void Solve2();

		/**
		* ������ ������� ������ ������ �����-�����
		*/
		virtual
			void Solve3();

		/**
		* ������ ��������� ������ ������ �����-�����
		*/
		virtual
			void Solve4();

		/**
		* ������ ����� ������ ������ �����-�����
		*/
		virtual
			void Solve5();

		virtual
			void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
#pragma endregion

		double df[12]; // debug
	};
}