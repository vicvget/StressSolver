#pragma once

#include "StressStrainCppIterativeSolver.h"

#include <vector>

using std::string;
using std::vector;

namespace Stress
{

/** Унаследованный код на ФОРТРАНе
*/
class StressStrainCppIterativeSolverAVX
	:
		public StressStrainCppIterativeSolver
{
public:

	// создает объект с заданными параметрами
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