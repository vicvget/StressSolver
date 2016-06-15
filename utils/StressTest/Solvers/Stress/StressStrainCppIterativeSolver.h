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

/** Унаследованный код на ФОРТРАНе
*/
class StressStrainCppIterativeSolver
	:
		public StressStrainCppSolver
{
public:

	// создает объект с заданными параметрами
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

	/** Получить смещения
	* @param data - массив для записи смещений как скалярного параметра
	*/
	virtual
	void GetDisplacement
		(
			float* data
		);

	/** Получить напряжения по первой теории прочности
	* @param data - массив для записи напряжений как скалярного параметра
	*/
	virtual
	void GetStressesByFirstTheoryOfStrength
		(
			float* data
		);

	/** Получить напряжения по von Mises
	* @param data - массив для записи напряжений как скалярного параметра
	*/
	virtual
	void GetStressesByVonMises
		(
			float* data
		);

#pragma endregion

	double df[12]; // debug

protected:
	// номер итерации расчетного цикла

	virtual void CalculateForces();

	int _iterationNumber;
	void ApplyBoundary();
	void ApplyMass();

};
}