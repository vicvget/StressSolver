//#ifdef USE_KNC
#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>
#include <immintrin.h>

using std::string;
using std::vector;

namespace Stress
{

	class StressStrainCppIterativeSolverKNC
		:
		public StressStrainCppIterativeSolver
	{
		//__m512d timeStep;
		//__m512d timeStep2;
		//__m512d timeStep4;
		//__m512d constantD2;
		//__m512d constantD6;
		const int regSize = 8;

	public:

		//virtual void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		// ������� ������ � ��������� �����������
		StressStrainCppIterativeSolverKNC
			(
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
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
//#endif