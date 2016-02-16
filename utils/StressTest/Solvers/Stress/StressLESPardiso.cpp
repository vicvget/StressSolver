#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_pardiso.h>
/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif

/** Обертка для расчетного модуля PARDISO
*/
class PardisoSolver
{
	MKL_INT mtype; /* Real symmetric matrix */
	/* RHS and solution vectors. */
	//double b[8], x[8];
	MKL_INT nrhs; /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */

	int* _ia; 
	int* _ja; 
	double* _b; 
	double* _a;
	double* _x;
	int _n;

	void Init(int* ia, int* ja, double* b, double* a, double* x, int n)
	{
		_n=n;
		_ia=ia;
		_ja=ja;
		_b=b;
		_a=a;
		_x=x;

		mtype = -2;
		nrhs = 1;
		for (i = 0; i < 64; i++) 
		{
			iparm[i] = 0;
		}
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
		maxfct = 1; /* Maximum number of numerical factorizations. */
		mnum = 1; /* Which factorization to use. */
		//msglvl = 1; /* Print statistical information in file */
		error = 0; /* Initialize error flag */
		/* -------------------------------------------------------------------- */
		/* .. Initialize the internal solver memory pointer. This is only */
		/* necessary for the FIRST call of the PARDISO solver. */
		/* -------------------------------------------------------------------- */
		for (i = 0; i < 64; i++) {
			pt[i] = 0;
		}
		/* -------------------------------------------------------------------- */
		/* .. Reordering and Symbolic Factorization. This step also allocates */
		/* all memory that is necessary for the factorization. */
		/* -------------------------------------------------------------------- */
		phase = 11;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&_n, _a, _ia, _ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			printf("\nERROR during symbolic factorization: %d", error);
			exit(1);
		}
		//printf("\nReordering completed ... ");
		//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
		//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
		/* -------------------------------------------------------------------- */
		/* .. Numerical factorization. */
		/* -------------------------------------------------------------------- */
		phase = 22;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&_n, _a, _ia, _ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			exit(2);
		}
	}



	bool SolveSLES() 
	{

		//	printf("\nFactorization completed ... ");
		/* -------------------------------------------------------------------- */
		/* .. Back substitution and iterative refinement. */
		/* -------------------------------------------------------------------- */
		phase = 33;
		iparm[7] = 2; /* Max numbers of iterative refinement steps. */
		/* Set right hand side to one. */
		//for (i = 0; i < n; i++) {
		//b[i] = 1;
		//}
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&_n, _a, _ia, _ja, &idum, &nrhs,
			iparm, &msglvl, _b, _x, &error);
		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			return false;
		}
		return true;
	}
	//printf("\nSolve completed ... ");
	//printf("\nThe solution of the system is: ");
	//for (i = 0; i < n; i++) {
	//	printf("\n x [%d] = % f", i, x[i] );
	//}
	//printf ("\n");
	void Dispose() 
	{
		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
		phase = -1; /* Release internal memory. */
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&_n, &ddum, _ia, _ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
	}

};