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
* Деструктор
*/
PardisoSolver::~PardisoSolver()
{
	if (_isAllocated)
	{
		Dispose();
	}
}

/**
* Инициализировать решатель MKL PARDISO
* @param n - количество уравнений в рассчитываемой системе
* @param ia - массив индексов первых ненулевых элементов в очередном ряду в массиве a
* @param ja - массив индексов столбцов для соответствующих элементов массива a
* @param a - массив ненулевых элементов матрицы A
* @param mtype - тип матрицы A
* @param isCSC - признак, передана ли матрица A в формате CSC (true) или CSR (false)
* @param useCStyleIndexing - признак, использовать ли индексацию в стиле C [с 0] (true) или нет (false)
* @return признак успешной (true) или неуспешной (false) инициализации
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
	// TODO: удалить после завершения отладки
	_iparm[26] = 1; // check whether column indices are sorted in increasing order within each row
	_iparm[34] = useCStyleIndexing ? 1 : 0; // C-style indexing

	// TODO: удалить старый код заполнения массива iparm всесторонней проверки нового
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
* Освобождение памяти, выделенной под нужды решателя MKL PARDISO
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
* Произвести расчет системы линейных уравнений решателем MKL PARDISO
* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
* @param x - массив, в который записываются результаты расчета
* @return признак успешного (true) или неуспешного (false) расчета
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
* Преобразовать матрицу из формата CSC (Compressed Sparse Column) в CSR (Compressed Sparse Row)
* @param n - размер преобразумой матрицы
* @param iaCSC - массив индексов первых ненулевых элементов в очередном ряду в массиве aCSC
* @param jaCSC - массив индексов столбцов для соответствующих элементов массива aCSC
* @param aCSC - массив ненулевых элементов матрицы в формате CSC
* @param iaCSR - массив индексов первых ненулевых элементов в очередном ряду в массиве aCSR
* @param jaCSR - массив индексов столбцов для соответствующих элементов массива aCSR
* @param aCSR - массив ненулевых элементов матрицы в формате CSR
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
* Вызвать метод pardiso с членами данного класса в качестве параметров
* @param phase - фаза работы решателя pardiso
* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
* @param x - массив, в который записываются результаты расчета
* @return 0 в случае успеха, код ошибки (< 0) в случае неудачи
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
* Произвести символьную и численную факторизацию матрицы A
* @return признак успешной (true) или неуспешной (false) факторизации
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
* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
* @param x - массив, в который записываются результаты расчета
* @return признак успеха (true) или неудачи (false)
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