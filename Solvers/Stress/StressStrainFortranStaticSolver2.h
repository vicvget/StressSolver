#ifndef STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H

#define STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H


#include "StressStrainFortranIterativeSolver.h"


/**
* Статический решатель напряженно-деформированного состояния тела на основе матрицы жесткости
*/
class StressStrainFortranStaticSolver2
	:
		public StressStrainFortranIterativeSolver
{
public:

	// создает объект с заданными параметрами
	StressStrainFortranStaticSolver2
		(
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
			double gridStep,
			double timeStep,
			int numThreads
		);

#pragma region overriden

	virtual
	void Solve
		(
			const int nIteratons
		);

	/**
	* Решить СЛАУ на основе матрицы жесткости
	* @param stiffnessMatrixCSC - матрица жесткости в формате CSC
	* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости
	* @return рассчитанный вектор результатов СЛАУ на основе матрицы жесткости
	*/
	static
	StiffnessRHSVector SolveSystemOfLinearEquationsForStiffness
		(
			const StiffnessMatrixPackedInCSCFormat& stiffnessMatrixCSC,
			const StiffnessRHSVector& stiffnessRHSVector
		);

#pragma endregion


private:

	/**
	* Сформировать СЛАУ на основе матрицы жесткости
	* @param stiffnessMatrixCSC - сформированная матрица жесткости в формате CSC
	* @param stiffnessRHSVector - сформированный вектор правых частей СЛАУ на основе матрицы жесткости
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void CreateSystemOfLinearEquationsForStiffness
		(
			StiffnessMatrixPackedInCSCFormat& stiffnessMatrixCSC,
			StiffnessRHSVector& stiffnessRHSVector,
			double stiffnessEpsilon = StiffnessEpsilon
		);

	/**
	* Использовать решение СЛАУ на основе матрицы жесткости для заполнения соответствующих массивов решателя
	* @param stiffnessResultsVector - вектор результатов СЛАУ на основе матрицы жесткости
	*/
	void UseSolutionOfSystemOfLinearEquationsForStiffness
		(
			const StiffnessRHSVector& stiffnessResultsVector
		);

};

#endif // STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H