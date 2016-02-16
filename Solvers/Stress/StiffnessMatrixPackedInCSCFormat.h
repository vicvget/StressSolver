#ifndef STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H

#define STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H


#include "StiffnessMatrixIO.h"

#include <vector>


/**
* Класс для хранения матрицы жесткости, упакованной в CSC (compressed sparse column) формате
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrixPackedInCSCFormat)
{
public:

	// Статические поля класса

	// суффикс, добавляемый к имени файла, в который записывается или из которого читается матрица жесткости
	static
	const std::string _fileNamePostfix;


	// Конструкторы и деструктор

	/**
	* Конструктор по умолчанию
	*/
	StiffnessMatrixPackedInCSCFormat() = default;

	/**
	* Конструктор
	* @param nodesCount количество узлов сеточного представления тела
	*/
	StiffnessMatrixPackedInCSCFormat
		(
			size_t nodesCount
		);


	// Селекторы

	/**
	* Получить размер полной (неупакованной) матрицы жесткости
	* @return размер полной (неупакованной) матрицы жесткости
	*/
	size_t GetMatrixSize() const;

	/**
	* Получить количество ненулевых элементов в матрице жесткости
	* @return количество ненулевых элементов в матрице жесткости
	*/
	size_t GetMatrixNonZeroElementsCount() const;

	/**
	* Получить указатель на массив элементов матрицы жесткости
	* @return указатель на массив элементов матрицы жесткости
	*/
	const double* GetElements() const;

	/**
	* Скопировать элементы матрицы жесткости в другой массив
	* @param elements - массив, в который производится копирвоание элементов матрицы жесткости
	*/
	void CopyElements
		(
			double* elements
		)	const;

	/**
	* Скопировать индексы строк элементов в матрице жесткости в другой массив
	* @param rowIds - массив, в который производится копирование индексов строк элементов в матрице жесткости
	*/
	void CopyRowIds
		(
			int* rowIds
		)	const;

	/**
	* Скопировать индексы первых элементов следующего столбца в предыдущем массиве в другой массив
	* @param columnsFirstIds - массив, в который производится копирвоание индексов первых элементов
	* следующего столбца в предыдущем массиве в другой массив
	*/
	void CopyColumnFirstElements
		(
			int* columnsFirstIds
		)	const;


	// Работа с элементами матрицы

	/**
	* Заполнить элемент матрицы жесткости
	* @param rowIndex индекс строки матрицы жесткости
	* @param columnIndex индекс столбца матрицы жесткости
	* @param elementValue значение заполняемого элемента матрицы жесткости
	*/
	void FillElement
		(
			size_t rowIndex,
			size_t columnIndex,
			double elementValue
		);

	/**
	* Заполнить элемент матрицы жесткости
	* @param rowIndex индекс ряда матрицы жесткости, состоящего из подматриц 6 x 6
	* @param rowDofIndex индекс строки подматрицы 6 x 6
	* @param columnIndex индекс колонки матрицы жесткости, состоящей из подматриц 6 x 6
	* @param columnDofIndex индекс столбца подматрицы 6 x 6
	* @param elementValue значение заполняемого элемента матрицы жесткости
	*/
	void FillElement
		(
			size_t rowIndex,
			size_t rowDofIndex,
			size_t columnIndex,
			size_t columnDofIndex,
			double elementValue
		);


	// Запись в файл / чтение из файла

	/**
	* Записать матрицу жесткости в текстовый файл
	* @param textFileName наименование текстового файла, в который производится запись матрицы жесткости
	*/
	void WriteToTextFile
		(
			const std::string& textFileName
		)	const;

	/**
	* Записать матрицу жесткости в текстовый файл формата GNU Octave
	* @param octaveTextFileName наименование текстового файла, в который производится запись матрицы жесткости
	*/
	void WriteToOctaveTextFile
		(
			const std::string& octaveTextFileName
		)	const;

	/**
	* Записать матрицу жесткости в двоичный файл
	* @param binaryFileName наименование двоичного файла, в который производится запись матрицы жесткости
	*/
	void WriteToBinaryFile
		(
			const std::string& binaryFileName
		)	const;

	/**
	* Считать матрицу жесткости из двоичного файла
	* @param binaryFileName наименование двоичного файла, из которого производится считывание матрицы жесткости
	*/
	void ReadFromBinaryFile
		(
			const std::string& binaryFileName
		);

private:

	// размер полной (неупакованной) матрицы жесткости
	size_t _matrixSize{};

	// значения элементов упакованной матрицы жесткости
	std::vector<double> _values;

	// индексы строк элементов в матрице жесткости
	std::vector<size_t> _rows;

	// индексы первых элементов следующего столбца в предыдущем массиве
	std::vector<size_t> _columnBeginnings = std::vector<size_t>{0};

};


// Явная конкретизация StiffnessMatrixIO

// для типа StiffnessMatrixPackedInCSCFormat
extern
template class StiffnessMatrixIO<StiffnessMatrixPackedInCSCFormat>;


#endif // STIFFNESS_MATRIX_PACKED_IN_CSC_FORMAT_H