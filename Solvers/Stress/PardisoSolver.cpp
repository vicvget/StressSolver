#ifdef USE_MKL

#ifdef _WIN64
#	pragma comment(lib, "mkl_intel_lp64.lib")
#else
#	pragma comment(lib, "mkl_intel_c.lib")
#endif

#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "mkl_core.lib")
#pragma comment(lib, "libiomp5md.lib")

#include "PardisoSolver.h"

#include <mkl_pardiso.h>
#include <mkl.h>
#include <iostream>


// class PardisoSolver

/**
* ����������
*/
PardisoSolver::~PardisoSolver()
{
	if (_isAllocated)
	{
		Dispose();
	}
}

/**
* ���������������� �������� MKL PARDISO
* @param n - ���������� ��������� � �������������� �������
* @param ia - ������ �������� ������ ��������� ��������� � ��������� ���� � ������� a
* @param ja - ������ �������� �������� ��� ��������������� ��������� ������� a
* @param a - ������ ��������� ��������� ������� A
* @param mtype - ��� ������� A
* @param isCSC - �������, �������� �� ������� A � ������� CSC (true) ��� CSR (false)
* @param useCStyleIndexing - �������, ������������ �� ���������� � ����� C [� 0] (true) ��� ��� (false)
* @return ������� �������� (true) ��� ���������� (false) �������������
*/
bool PardisoSolver::Init
	(
		MKL_INT n,
		MKL_INT* ia,
		MKL_INT* ja,
		double* a,
		MKL_INT mtype,
		bool isCSC,
		bool useCStyleIndexing
	)
{
	_n = n;
	_ia = ia;
	_ja = ja;
	_a = a;

	_mtype = mtype; // Matrix type

	pardisoinit(_pt, &mtype, _iparm);

	_iparm[11] = isCSC ? 2 : 0; // solve a transposed system A^T * X = B based on the factorization of the matrix A
	// TODO: ������� ����� ���������� �������
	_iparm[26] = 1; // check whether column indices are sorted in increasing order within each row
	_iparm[34] = useCStyleIndexing ? 1 : 0; // C-style indexing

	// TODO: ������� ������ ��� ���������� ������� iparm ������������ �������� ������
#if 0
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */

	iparm[27] = 1; /*matrix checker*/
	iparm[35] = 1; /*C-style indexing*/
#endif

	return true;
}

/**
* ������������ ������, ���������� ��� ����� �������� MKL PARDISO
*/
void PardisoSolver::Dispose()
{
	if (MKL_INT error = CallPardiso(-1))
	{
		std::cout << "ERROR during releasing memory: " << error << std::endl;

		return;
	}
	_isAllocated = false;
}

/**
* ���������� ������ ������� �������� ��������� ��������� MKL PARDISO
* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
* @param x - ������, � ������� ������������ ���������� �������
* @return ������� ��������� (true) ��� ����������� (false) �������
*/
bool PardisoSolver::Solve
	(
		double* b,
		double* x
	)
{
	if (Factorization())
	{
		return SolveBackward(b, x);
	}
	else
	{
		return false;
	}
}

/**
* ������������� ������� �� ������� CSC (Compressed Sparse Column) � CSR (Compressed Sparse Row)
* @param n - ������ ������������ �������
* @param iaCSC - ������ �������� ������ ��������� ��������� � ��������� ���� � ������� aCSC
* @param jaCSC - ������ �������� �������� ��� ��������������� ��������� ������� aCSC
* @param aCSC - ������ ��������� ��������� ������� � ������� CSC
* @param iaCSR - ������ �������� ������ ��������� ��������� � ��������� ���� � ������� aCSR
* @param jaCSR - ������ �������� �������� ��� ��������������� ��������� ������� aCSR
* @param aCSR - ������ ��������� ��������� ������� � ������� CSR
*/
// static
void PardisoSolver::ConvertFromCSC2CSR
	(
		MKL_INT n,
		MKL_INT* iaCSC,
		MKL_INT* jaCSC,
		double* aCSC,
		MKL_INT* iaCSR,
		MKL_INT* jaCSR,
		double* aCSR
	)
{
	int job[6]{1}; // convert from CSC to CSR
	int info;

	job[5] = 1; // fill 'elements' array
	mkl_dcsrcsc(job, &n, aCSR, jaCSR, iaCSR, aCSC, jaCSC, iaCSC, &info);
}

/**
* ������� ����� pardiso � ������� ������� ������ � �������� ����������
* @param phase - ���� ������ �������� pardiso
* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
* @param x - ������, � ������� ������������ ���������� �������
* @return 0 � ������ ������, ��� ������ (< 0) � ������ �������
*/
MKL_INT PardisoSolver::CallPardiso
	(
		MKL_INT phase,
		double* b,
		double* x
	)
{
	MKL_INT error{};

	PARDISO
		(
			_pt,
			&_maxfct,
			&_mnum,
			&_mtype,
			&phase,
			&_n,
			_a,
			_ia,
			_ja,
			nullptr,
			&_nrhs,
			_iparm,
			&_msglvl,
			b,
			x,
			&error
		);

	return error;
}

/**
* ���������� ���������� � ��������� ������������ ������� A
* @return ������� �������� (true) ��� ���������� (false) ������������
*/
bool PardisoSolver::Factorization()
{
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	if (MKL_INT error = CallPardiso(11))
	{
		std::cout << "ERROR during symbolic factorization: " << error << std::endl;

		return false;
	}
	_isAllocated = true;

	//printf("\nReordering completed ... ");
	//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	if (MKL_INT error = CallPardiso(22))
	{
		std::cout << "ERROR during numerical factorization: " << error << std::endl;

		return false;
	}

	return true;
}

/**
* Back substitution and iterative refinement
* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
* @param x - ������, � ������� ������������ ���������� �������
* @return ������� ������ (true) ��� ������� (false)
*/
bool PardisoSolver::SolveBackward
	(
		double* b,
		double* x
	)
{
	_iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	if (MKL_INT error = CallPardiso(33, b, x))
	{
		std::cout << "ERROR during solution: " << error << std::endl;

		return false;
	}

	return true;
}

#endif