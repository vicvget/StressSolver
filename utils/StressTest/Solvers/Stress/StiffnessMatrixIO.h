#ifndef STIFFNESS_MATRIX_IO_H

#define STIFFNESS_MATRIX_IO_H


#include "StiffnessMatrixFileFormat.h"

#include <string>


/**
* Класс для организации ввода/вывода матрицы жесткости в/из файлы/-ов различных форматов
*/
template
	<
		class Matrix // тип матрицы жесткости
	>
class StiffnessMatrixIO
{
public:

	// Запись в файл / чтение из файла

	/**
	* Записать матрицу жесткости в файл
	* @param outputFileName наименование файла, в который производится запись матрицы жесткости
	* @param fileFormats набор форматов файлов, в которые должна производиться запись
	*/
	void WriteToFile
		(
			const std::string& outputFileName,
			StiffnessMatrixFileFormats fileFormats
		)	const;


private:

	// внутренние типы

	// текущий тип класса для организации ввода/вывода матрицы жесткости
	using ThisType = StiffnessMatrixIO<Matrix>;

	// тип функции-члена класса матрицы жесткости, предназначенного для ее вывода
	typedef
	void (Matrix::*WriteMember)
		(
			const std::string& outputFileName // наименование файла, в который выводится матрица жесткости
		)	const;


	// статические переменные-члены

	// массив функций-членов класса матрицы жесткости, предназначенных для ее вывода
	static
	const WriteMember _writeMembers[];

	// массив строк, содержащий расширения файлов матриц жесткости
	static
	const std::string _fileFormatsExtensions[];

};


// объявить класс матрицы жесткости
#define DECLARE_STIFFNESS_MATRIX(Matrix) \
 \
class Matrix final \
	: \
		public StiffnessMatrixIO<Matrix>


#include "StiffnessMatrixIO.hpp"


#endif // STIFFNESS_MATRIX_IO_H