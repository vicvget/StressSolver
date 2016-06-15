#pragma once

#include "StressStrainCppSolver.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>


using std::string;
using std::vector;

namespace Stress
{

/** �������������� ��� �� ��������
*/
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
			double *nodes,
			int nNodes,
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
	virtual void SolveFull(const int nIteratons);;

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
	virtual
	void GetStressesByVonMises
		(
			float* data
		);

#pragma endregion

	double df[12]; // debug

protected:
	// ����� �������� ���������� �����

	virtual void CalculateForces();

	int _iterationNumber;
	void ApplyBoundary();
	void ApplyMass();

};
}