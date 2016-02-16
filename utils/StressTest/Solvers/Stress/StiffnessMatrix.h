#ifndef STIFFNESS_MATRIX_H

#define STIFFNESS_MATRIX_H


#include "AuxiliaryStressStuff.h"
#include "StiffnessMatrixIO.h"


#include <memory>
#include <string>


// матрица 6 x 6, которая используется для хранения элементов матрицы жесткости
using Matrix6x6 = double[FreedomsCount][FreedomsCount];

// указатель на матрицу 6 x 6
using PMatrix6x6 = std::unique_ptr<Matrix6x6[]>;

// двойной указатель на матрицу 6 x 6
using PPMatrix6x6 = std::unique_ptr<PMatrix6x6[]>;


/**
* Матрица жесткости, составленная из подматриц размера 6 x 6
* (количество таких подматриц равно _nodesCount x _nodesCount)
*/
DECLARE_STIFFNESS_MATRIX(StiffnessMatrix)
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
	StiffnessMatrix
		(
			size_t nodesCount
		);

	/**
	* Конструктор перемещения
	* @param otherStiffnessMatrix - матрица жесткости, из которой производится перемещение
	*/
	StiffnessMatrix
		(
			StiffnessMatrix&& otherStiffnessMatrix
		);


	// Работа с элементами матрицы

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


	// Доступ к элементам матрицы жесткости
	
	/**
	* Получить указатель на ряд матрицы жесткости с индексом rowIndex (ряд состоит из подматриц 6 x 6)
	* @param rowIndex индекс ряда матрицы жесткости
	* @return указатель на ряд матрицы жесткости с индексом rowIndex
	*/
	PMatrix6x6& operator[]
		(
			size_t rowIndex
		);


private:

	// количество узлов в сеточном представлении тела
	size_t _nodesCount{};

	// двойной указатель на матрицу 6 x 6, по которому хранится матрица жесткости
	PPMatrix6x6 _matrix{};


	// Вспомогательные функции

	/**
	* Выделить память для хранения матрицы жесткости
	* (по известному числу узлов в сеточном представлении _nodesCount)
	*/
	void Allocate();

};


// Явная конкретизация StiffnessMatrixIO

// для типа StiffnessMatrix
extern
template class StiffnessMatrixIO<StiffnessMatrix>;


#endif // STIFFNESS_MATRIX_H