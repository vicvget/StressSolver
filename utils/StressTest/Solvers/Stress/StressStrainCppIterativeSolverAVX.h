#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>

using std::string;
using std::vector;

namespace Stress
{

/** �������������� ��� �� ��������
*/
class StressStrainCppIterativeSolverAVX
	:
		public StressStrainCppIterativeSolver
{
public:

	// ������� ������ � ��������� �����������
	StressStrainCppIterativeSolverAVX
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
	~StressStrainCppIterativeSolverAVX();

#pragma region overriden

	virtual	void InitialSolve();

	virtual void Solve(const int nIteratons);
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

#pragma endregion

	double df[12]; // debug
protected:

	virtual	void CalculateForces();
};
}