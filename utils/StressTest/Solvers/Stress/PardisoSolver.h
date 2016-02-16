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
* Обертка для расчетного модуля MKL PARDISO.
* Производит расчет системы линейных уравнений вида A * X = B
*/
class PardisoSolver final
{

public:

	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	PardisoSolver() = default;

	/**
	* Конструктор копирования
	*/
	PardisoSolver(const PardisoSolver&) = delete;

	/**
	* Деструктор
	*/
	~PardisoSolver();


	// Инициализация и деинициализация

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
	* Освобождение памяти, выделенной под нужды решателя MKL PARDISO
	*/
	void Dispose();


	// Расчет

	/**
	* Произвести расчет системы линейных уравнений решателем MKL PARDISO
	* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
	* @param x - массив, в который записываются результаты расчета
	* @return признак успешного (true) или неуспешного (false) расчета
	*/
	bool Solve
		(
			double* b,
			double* x
		);


	// Вспомогательные статические функции

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

	// признак, выделена ли память под внутренние нужды решателя MKL PARDISO
	bool _isAllocated{};

	// MKL PARDISO control parameters

	// массив параметров решателя MKL PARDISO
	MKL_INT _iparm[64]{};

	// maximum number of factors with identical sparsity structure that must be kept in memory at the same time
	MKL_INT _maxfct{1};

	// define which matrix to factorize
	MKL_INT _mnum{1};

	// message level information
	MKL_INT _msglvl{}; // don't print statistical information in file

	// параметры рассчитываемой системы уравнений

	// количество уравнений в рассчитываемой системе
	MKL_INT _n;

	// массив индексов первых ненулевых элементов в очередном ряду в массиве a
	MKL_INT* _ia;

	// массив индексов столбцов для соответствующих элементов массива a
	MKL_INT* _ja;

	// массив ненулевых элементов матрицы A
	double* _a;


	// Вспомогательные функции

	/**
	* Вызвать метод pardiso с членами данного класса в качестве параметров
	* @param phase - фаза работы решателя pardiso
	* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
	* @param x - массив, в который записываются результаты расчета
	* @return 0 в случае успеха, код ошибки (< 0) в случае неудачи
	*/
	MKL_INT CallPardiso
		(
			MKL_INT phase,
			double* b = nullptr,
			double* x = nullptr
		);

	/**
	* Произвести символьную и численную факторизацию матрицы A
	* @return признак успешной (true) или неуспешной (false) факторизации
	*/
	bool Factorization();

	/**
	* Back substitution and iterative refinement
	* @param b - массив, в котором находятся значения вектора правых частей системы уравнений
	* @param x - массив, в который записываются результаты расчета
	* @return признак успеха (true) или неудачи (false)
	*/
	bool SolveBackward
		(
			double* b,
			double* x
		);

};

#endif

#endif // PARDISO_SOLVER_H