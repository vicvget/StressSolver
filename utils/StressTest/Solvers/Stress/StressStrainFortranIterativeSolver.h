#ifndef StressStrainFortranIterativeSolverH

#define StressStrainFortranIterativeSolverH


#include "StressStrainFortranSolver.h"
#include "StiffnessMatrices.h"
#include "StiffnessRHS.h"

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>


// наименование файла, в который выводится матрица жесткости, по умолчанию
const std::string StiffnessMatrixFileName = "StiffnessMatrix";

// наименование файла, в который выводится вектор правых частей СЛАУ на основе матрицы жесткости, по умолчанию
const std::string StiffnessRHSFileName = "StiffnessRHS";

// число, использующееся в численном нахождении производной при формировании матрицы жесткости
const double StiffnessEpsilon = 1e-5;


/**
* Вспомогательная структура, использующаяся для идентификации типа
*/
template
	<
		typename Type // идентифицируемый тип
	>
struct Identity
{
};


/** Унаследованный код на ФОРТРАНе
*/
class StressStrainFortranIterativeSolver
	:public StressStrainFortranSolver
{
public:

	// создает объект с заданными параметрами
	StressStrainFortranIterativeSolver
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

#pragma region overriden

	virtual
	void Solve
		(
			const int nIteratons
		);

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


//protected:

	// используемый тип матрицы жесткости (или кортеж матриц жетскости)
	using StiffnessMatrixType = StiffnessMatrixPackedInCSCFormat;


	// номер итерации расчетного цикла
	int _iterationNumber;


	void pravsubfl();

	void FindStressStrainMatrix
		(
			double (&matrix)[6][6]
		);

	void ApplyBoundary();

	void ApplyMass();

	/**
	* Сформировать СЛАУ на основе матрицы жесткости
	* @param stiffnessMatrixFileName - общее наименование файлов для записи формируемой матрицы жесткости
	* в ряде форматов
	* @param stiffnessMatrixFileFormats - список форматов файлов, в которые будет записана
	* формируемая матрица жесткости
	* @param stiffnessRHSFileName - общее наименование файлов для записи формируемого вектора
	* правых частей СЛАУ на основе матрицы жесткости
	* @param writeStiffnessRHSVectorToTextFileFlag - флаг, следует ли записывать формируемый вектор
	* правых частей СЛАУ на основе матрицы жесткости в текстовый файл
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void CreateSystemOfLinearEquationsForStiffness
		(
			const string& stiffnessMatrixFileName = StiffnessMatrixFileName,
			StiffnessMatrixFileFormats stiffnessMatrixFileFormats = //TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat::BinaryFile),
				TO_INTEGRAL_TYPE(StiffnessMatrixFileFormat::BinaryFile | StiffnessMatrixFileFormat::TextFile),
			const string& stiffnessRHSFileName = StiffnessRHSFileName,
			//bool writeStiffnessRHSVectorToTextFileFlag = false,
			bool writeStiffnessRHSVectorToTextFileFlag = true,
			double stiffnessEpsilon = StiffnessEpsilon
		);

	/**
	* Сформировать список наборов флагов закрепленных степеней свободы элементов сеточного представления
	* @return формируемый список наборов флагов закрепленных степеней свободы элементов сеточного представления
	*/
	IsSealedFlagsList FormIsSealedFlagsList() const;

	/**
	* Сформировать матрицу жесткости
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessMatrixFileName - общее наименование файлов для записи формируемой матрицы жесткости
	* в ряде форматов
	* @param stiffnessMatrixFileFormats - список форматов файлов, в которые будет записана
	* формируемая матрица жесткости
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void CreateStiffnessMatrix
		(
			const IsSealedFlagsList& isSealedFlagsList,
			const string& stiffnessMatrixFileName,
			StiffnessMatrixFileFormats stiffnessMatrixFileFormats,
			double stiffnessEpsilon
		);

	/**
	* Сформировать матрицу жесткости
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	* @return сформированная матрица жесткости
	*/
	StiffnessMatrixType CreateStiffnessMatrix
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* Сформировать матрицу жесткости, упакованную в CSC (compressed sparse column) формате
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	* @param вспомогательный параметр, необходимый для выведения типа UsedStiffnessMatrixType
	* @return сформированная матрица жесткости, упакованная в CSC (compressed sparse column) формате
	*/
	template
		<
			class UsedStiffnessMatrixType = StiffnessMatrixType // используемый тип матрицы жесткости
		>
	std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrixPackedInCSCFormat>::value, StiffnessMatrixPackedInCSCFormat>
	CreateStiffnessMatrixPackedInCSCFormat
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon,
			Identity<UsedStiffnessMatrixType> = {}
		);

	/**
	* Сформировать матрицу жесткости, упакованную в CSC (compressed sparse column) формате
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	* @param вспомогательный параметр, необходимый для выведения типа UsedStiffnessMatrixType
	* @return сформированная матрица жесткости, упакованная в CSC (compressed sparse column) формате
	*/
	template
		<
			class UsedStiffnessMatrixType = StiffnessMatrixType // используемый тип матрицы жесткости
		>
	std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrices>::value, StiffnessMatrixPackedInCSCFormat>
	CreateStiffnessMatrixPackedInCSCFormat
		(
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon,
			Identity<UsedStiffnessMatrixType> = {}
		);

	/**
	* Заполнить матрицу жесткости
	* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void FindStiffnessMatrix
		(
			StiffnessMatrixType& stiffnessMatrix,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* Заполнить матрицу жесткости значениями всех соседних элементов к элементу с данным индексом
	* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
	* @param nodeId - индекс данного элемента
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void FindStiffnessMatrixForAdjElements
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* Заполнить матрицу жесткости значениями соседнего элемента к элементу с данным индексом
	* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
	* @param nodeId - индекс данного элемента
	* @param adjNodeId - индекс соседнего элемента
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void FindStiffnessMatrixForElement
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			int adjNodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* Заполнить матрицу жесткости значениями элемента, некоторые степени свободы которого закреплены
	* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
	* @param nodeId - индекс элемента, некоторые степени свободы которого закреплены
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	*/
	void FindStiffnessMatrixForSealedElement
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			const IsSealedFlagsList& isSealedFlagsList
		);

	/**
	* Заполнить матрицу жесткости значениями соседнего элемента к элементу с данным индексом
	* и данной степенью свободы
	* @param stiffnessMatrix - матрица жесткости, которая заполняется значениями
	* @param nodeId - индекс данного элемента
	* @param dof - индекс данной степени свободы
	* @param adjNodeId - индекс соседнего элемента
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void FindStiffnessMatrixForDof
		(
			StiffnessMatrixType& stiffnessMatrix,
			int nodeId,
			int dof,
			int adjNodeId,
			const IsSealedFlagsList& isSealedFlagsList,
			double stiffnessEpsilon
		);

	/**
	* Обновить значения массива ускорений для элемента с данным индексом, изменив значение приложенной
	* по направлению changedDof силы к элементу с индексом changedNodeId
	* @param changedNodeId - индекс элемента, к которому прикладывается сила
	* @param changedDof - индекс степени свободы, по которой прикладывается сила
	* @param nodeId - индекс данного элемента
	* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
	* при формировании матрицы жесткости
	*/
	void UpdateForcesForNode
		(
			int changedNodeId,
			int changedDof,
			int nodeId,
			double stiffnessEpsilon
		);

	/**
	* Сформировать вектор правых частей СЛАУ на основе матрицы жесткости
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @param stiffnessRHSFileName - общее наименование файлов для записи формируемого вектора
	* правых частей СЛАУ на основе матрицы жесткости
	* @param writeStiffnessRHSVectorToTextFileFlag - флаг, следует ли записывать формируемый вектор
	* правых частей СЛАУ на основе матрицы жесткости в текстовый файл
	*/
	void CreateStiffnessRHSVector
		(
			const IsSealedFlagsList& isSealedFlagsList,
			const string& stiffnessRHSFileName,
			bool writeStiffnessRHSVectorToTextFileFlag
		)	const;

	/**
	* Сформировать вектор правых частей СЛАУ на основе матрицы жесткости
	* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
	* элементов сеточного представления
	* @return сформированный вектор правых частей СЛАУ на основе матрицы жесткости, заполняемый значениями
	*/
	StiffnessRHSVector CreateStiffnessRHSVector
		(
			const IsSealedFlagsList& isSealedFlagsList
		)	const;

	/**
	* Заполнить значениями вектор правых частей СЛАУ на основе матрицы жесткости
	* @param stiffnessRHSVector - вектор правых частей СЛАУ на основе матрицы жесткости, заполняемый значениями
	*/
	void FindStiffnessRHSVector
		(
			StiffnessRHSVector& stiffnessRHSVector
		)	const;

};

/**
* Сформировать матрицу жесткости, упакованную в CSC (compressed sparse column) формате
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
* @param вспомогательный параметр, необходимый для выведения типа UsedStiffnessMatrixType
* @return сформированная матрица жесткости, упакованная в CSC (compressed sparse column) формате
*/
template
	<
		class UsedStiffnessMatrixType
	>
std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrixPackedInCSCFormat>::value, StiffnessMatrixPackedInCSCFormat>
StressStrainFortranIterativeSolver::CreateStiffnessMatrixPackedInCSCFormat
	(
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon,
		Identity<UsedStiffnessMatrixType>
	)
{
	return CreateStiffnessMatrix(isSealedFlagsList, stiffnessEpsilon);
}

/**
* Сформировать матрицу жесткости, упакованную в CSC (compressed sparse column) формате
* @param isSealedFlagsList - список наборов флагов закрепленных степеней свободы
* элементов сеточного представления
* @param stiffnessEpsilon - число, использующееся в численном нахождении производной
* при формировании матрицы жесткости
* @param вспомогательный параметр, необходимый для выведения типа UsedStiffnessMatrixType
* @return сформированная матрица жесткости, упакованная в CSC (compressed sparse column) формате
*/
template
	<
		class UsedStiffnessMatrixType
	>
std::enable_if_t<std::is_same<UsedStiffnessMatrixType, StiffnessMatrices>::value, StiffnessMatrixPackedInCSCFormat>
StressStrainFortranIterativeSolver::CreateStiffnessMatrixPackedInCSCFormat
	(
		const IsSealedFlagsList& isSealedFlagsList,
		double stiffnessEpsilon,
		Identity<UsedStiffnessMatrixType>
	)
{
	return std::get<0>(CreateStiffnessMatrix(isSealedFlagsList, stiffnessEpsilon).GetTuple());
}

#endif