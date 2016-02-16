#ifndef PARDISO_SOLVER_H

#define PARDISO_SOLVER_H

#ifdef USE_MKL

#if defined(_WIN32) || defined(_WIN64)
#	define pardiso_ PARDISO
#else
#	define PARDISO pardiso_
#endif

#if defined(MKL_ILP64)
#	define MKL_INT long long
#else
#	define MKL_INT int
#endif


/**
* ������� ��� ���������� ������ MKL PARDISO.
* ���������� ������ ������� �������� ��������� ���� A * X = B
*/
class PardisoSolver final
{

public:

	// ������������ � ����������

	/**
	* ����������� �� ���������
	*/
	PardisoSolver() = default;

	/**
	* ����������� �����������
	*/
	PardisoSolver(const PardisoSolver&) = delete;

	/**
	* ����������
	*/
	~PardisoSolver();


	// ������������� � ���������������

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
	bool Init
		(
			MKL_INT n,
			MKL_INT* ia,
			MKL_INT* ja,
			double* a,
			MKL_INT mtype = 11, // real and nonsymmetric matrix
			bool isCSC = true, // matrix in CSC format
			bool useCStyleIndexing = true // use C-style indexing
		);

	/**
	* ������������ ������, ���������� ��� ����� �������� MKL PARDISO
	*/
	void Dispose();


	// ������

	/**
	* ���������� ������ ������� �������� ��������� ��������� MKL PARDISO
	* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
	* @param x - ������, � ������� ������������ ���������� �������
	* @return ������� ��������� (true) ��� ����������� (false) �������
	*/
	bool Solve
		(
			double* b,
			double* x
		);


	// ��������������� ����������� �������

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
	static
	void ConvertFromCSC2CSR
		(
			MKL_INT n,
			MKL_INT* iaCSC,
			MKL_INT* jaCSC,
			double* aCSC,
			MKL_INT* iaCSR,
			MKL_INT* jaCSR,
			double* aCSR
		);


private:

	// matrix type 
	MKL_INT _mtype;

	// number of right hand sides
	MKL_INT _nrhs{1};

	// internal solver memory pointer pt
	void* _pt[64]{};

	// �������, �������� �� ������ ��� ���������� ����� �������� MKL PARDISO
	bool _isAllocated{};

	// MKL PARDISO control parameters

	// ������ ���������� �������� MKL PARDISO
	MKL_INT _iparm[64]{};

	// maximum number of factors with identical sparsity structure that must be kept in memory at the same time
	MKL_INT _maxfct{1};

	// define which matrix to factorize
	MKL_INT _mnum{1};

	// message level information
	MKL_INT _msglvl{}; // don't print statistical information in file

	// ��������� �������������� ������� ���������

	// ���������� ��������� � �������������� �������
	MKL_INT _n;

	// ������ �������� ������ ��������� ��������� � ��������� ���� � ������� a
	MKL_INT* _ia;

	// ������ �������� �������� ��� ��������������� ��������� ������� a
	MKL_INT* _ja;

	// ������ ��������� ��������� ������� A
	double* _a;


	// ��������������� �������

	/**
	* ������� ����� pardiso � ������� ������� ������ � �������� ����������
	* @param phase - ���� ������ �������� pardiso
	* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
	* @param x - ������, � ������� ������������ ���������� �������
	* @return 0 � ������ ������, ��� ������ (< 0) � ������ �������
	*/
	MKL_INT CallPardiso
		(
			MKL_INT phase,
			double* b = nullptr,
			double* x = nullptr
		);

	/**
	* ���������� ���������� � ��������� ������������ ������� A
	* @return ������� �������� (true) ��� ���������� (false) ������������
	*/
	bool Factorization();

	/**
	* Back substitution and iterative refinement
	* @param b - ������, � ������� ��������� �������� ������� ������ ������ ������� ���������
	* @param x - ������, � ������� ������������ ���������� �������
	* @return ������� ������ (true) ��� ������� (false)
	*/
	bool SolveBackward
		(
			double* b,
			double* x
		);

};

#endif

#endif // PARDISO_SOLVER_H