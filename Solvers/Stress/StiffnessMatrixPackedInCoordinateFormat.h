#ifndef STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H

#define STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H


#include "StiffnessMatrixIO.h"

#include <vector>
#include <string>


/**
* Элемент матрицы жесткости, упакованной в координатном формате
*/
struct MatrixElement final
{

	// индекс строки матрицы жесткости, в которой расположен данный элемент
	int rowIndex;

	// индекс столбца матрицы жесткости, в котором расположен данный элемент
	int columnIndex;

	// значение элемента матрицы жесткости с данными индексами
	double elementValue;

};


/**
* Класс для хранения матрицы жесткости, упакованной в координатном формате
* (для каждого элемента полной матрицы хранятся его индексы и значение)
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrixPackedInCoordinateFormat)
{
public:

	// Статические поля класса

	// суффикс, добавляемый к имени файла, в который записывается или из которого читается матрица жесткости
	static
	const std::string _fileNamePostfix;


	// Конструкторы и деструктор

	/**
	* Конструктор
	* @param nodesCount количество узлов сеточного представления тела
	*/
	StiffnessMatrixPackedInCoordinateFormat
		(
			size_t nodesCount
		);


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

	// элементы упакованной матрицы жесткости
	std::vector<MatrixElement> _matrixElements;

};


// Явная конкретизация StiffnessMatrixIO

// для типа StiffnessMatrixPackedInCoordinateFormat
extern
template class StiffnessMatrixIO<StiffnessMatrixPackedInCoordinateFormat>;


#endif // STIFFNESS_MATRIX_PACKED_IN_COORDINATE_FORMAT_H