#pragma once

#include "StressStrainCppSolver.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <immintrin.h>


using std::string;
using std::vector;

namespace Stress
{

	class StressStrainCppIterativeSolver
		:
		public StressStrainCppSolver
	{
	public:

		// ������� ������ � ��������� �����������
		StressStrainCppIterativeSolver
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
			~StressStrainCppIterativeSolver();

#pragma region overriden

		virtual
			void InitialSolve
			(
			);

		virtual	void Solve(const int nIteratons);
		virtual void SolveFull(const int nIteratons);

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

		/** �������� ��������
		* @param data - ������ ��� ������ �������� ��� ���������� ���������
		*/
		virtual
			void GetDisplacement
			(
			float* data
			);

		/** �������� ���������� �� ������ ������ ���������
		* @param data - ������ ��� ������ ���������� ��� ���������� ���������
		*/
		virtual
			void GetStressesByFirstTheoryOfStrength
			(
			float* data
			);

		/** �������� ���������� �� von Mises
		* @param data - ������ ��� ������ ���������� ��� ���������� ���������
		*/
		virtual	void GetStressesByVonMises(float* data);
		virtual
			void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		virtual
			void CalculateStrainsUa(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		/** �������� ���������� �� von Mises
				* @param data - ������ ��� ������ ���������� ��� ���������� ���������
				*/
		virtual void GetStressesX(float* data);
		virtual void GetStressesY(float* data);
		virtual void GetStressesZ(float* data);

		virtual void GetStressesXY(float* data);
		virtual void GetStressesXZ(float* data);
		virtual void GetStressesYZ(float* data);

#pragma endregion

		double df[12]; // debug

	protected:
		// ����� �������� ���������� �����
		int _iterationNumber;

	};
}