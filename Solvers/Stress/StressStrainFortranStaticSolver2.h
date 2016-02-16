#ifndef STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H

#define STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H


#include "StressStrainFortranIterativeSolver.h"


/**
* ����������� �������� ����������-���������������� ��������� ���� �� ������ ������� ���������
*/
class StressStrainFortranStaticSolver2
	:
		public StressStrainFortranIterativeSolver
{
public:

	// ������� ������ � ��������� �����������
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
	* ������ ���� �� ������ ������� ���������
	* @param stiffnessMatrixCSC - ������� ��������� � ������� CSC
	* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
	* @return ������������ ������ ����������� ���� �� ������ ������� ���������
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
	* ������������ ���� �� ������ ������� ���������
	* @param stiffnessMatrixCSC - �������������� ������� ��������� � ������� CSC
	* @param stiffnessRHSVector - �������������� ������ ������ ������ ���� �� ������ ������� ���������
	* @param stiffnessEpsilon - �����, �������������� � ��������� ���������� �����������
	* ��� ������������ ������� ���������
	*/
	void CreateSystemOfLinearEquationsForStiffness
		(
			StiffnessMatrixPackedInCSCFormat& stiffnessMatrixCSC,
			StiffnessRHSVector& stiffnessRHSVector,
			double stiffnessEpsilon = StiffnessEpsilon
		);

	/**
	* ������������ ������� ���� �� ������ ������� ��������� ��� ���������� ��������������� �������� ��������
	* @param stiffnessResultsVector - ������ ����������� ���� �� ������ ������� ���������
	*/
	void UseSolutionOfSystemOfLinearEquationsForStiffness
		(
			const StiffnessRHSVector& stiffnessResultsVector
		);

};

#endif // STRESS_STRAIN_FORTRAN_STATIC_SOLVER2_H