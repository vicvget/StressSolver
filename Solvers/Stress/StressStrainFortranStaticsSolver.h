#ifndef StressStrainFortranStaticsSolverH

#define StressStrainFortranStaticsSolverH

#include "StressStrainFortranSolver.h"
#include "PardisoSolver.h"
#include "CsrSymmetricMatrix.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>


using std::string;
using std::vector;


/** �������������� ��� �� ��������
*/
class StressStrainFortranStaticsSolver
	:
		public StressStrainFortranSolver
{

public:

	// ������� ������ � ��������� �����������
	StressStrainFortranStaticsSolver
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

	virtual
	~StressStrainFortranStaticsSolver();
	
#pragma region overriden

	virtual
	void Solve
		(
			const int nIteratons
		);

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


private:
#ifdef USE_MKL
	PardisoSolver _pardisoSolver;
#endif
	CsrSymmetricMatrix _mtx;

	void ApplyBoundary();
	bool SolveLes();	
	void InitSolveLES();

	void MakeKmtx1d(int n, int* &ia, int* &ja, double* &a);
	void MakeKmtx1d(int n,CsrSymmetricMatrix& mtx);
	//void MakeKmtx3d(int n, int* &ia, int* &ja, double* &a);
	//void MakeKmtx3d(int n, int n2,
	//		int* &ia, int* &ja, double* &ak, 
	//		double a[6][6], double c[6][6]);
	void Genmatl(int n1, int n2, double* a, double* c);	
	void pravsubfl();

	void ApplyForceBoundary(const BoundaryParams& boundaryParams);

	void ApplySealedBoundary(const BoundaryParams& boundaryParams);

};
#endif