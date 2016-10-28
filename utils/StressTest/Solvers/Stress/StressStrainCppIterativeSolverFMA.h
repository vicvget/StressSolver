#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>
#include <immintrin.h>

using std::string;
using std::vector;

namespace Stress
{

class StressStrainCppIterativeSolverFMA
	:
		public StressStrainCppIterativeSolver
{
public:
	const int regSize = 4;

	// создает объект с заданными параметрами
	StressStrainCppIterativeSolverFMA
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
	~StressStrainCppIterativeSolverFMA();

#pragma region overriden

	virtual	void InitialSolve();

	virtual void Solve(const int nIteratons);
	virtual	void SolveFull(const int nIteratons);

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

#pragma endregion

	double df[12]; // debug
protected:

	virtual	void CalculateForces();
};
}