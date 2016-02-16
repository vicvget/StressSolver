#ifndef STIFFNESS_MATRIX_IO_HPP

#define STIFFNESS_MATRIX_IO_HPP


#include "StiffnessMatrixIO.h"


/**
* Записать матрицу жесткости в файл
* @param outputFileName наименование файла, в который производится запись матрицы жесткости
* @param fileFormats набор форматов файлов, в которые должна производиться запись
*/
template
	<
		class Matrix // тип матрицы жесткости
	>
void StiffnessMatrixIO<Matrix>::WriteToFile
	(
		const std::string& outputFileName,
		StiffnessMatrixFileFormats fileFormats
	)	const
{
	const Matrix& thisMatrix = static_cast<const Matrix&>(*this);
	const std::string fullOutputFileName = outputFileName + Matrix::_fileNamePostfix;
	StiffnessMatrixFileFormats currentFormat = 1;

	for
		(
			size_t formatIndex = 0;
			formatIndex < StiffnessMatrixFileFormatsCount;
			++formatIndex,
			currentFormat <<= 1
		)
	{
		const std::string outputFileNameWithExtension =
			fullOutputFileName + '.' + _fileFormatsExtensions[formatIndex];

		if (currentFormat & fileFormats)
		{
			(thisMatrix.*_writeMembers[formatIndex])(outputFileNameWithExtension);
		}
	}
}


// формат файла матрицы жесткости в качестве элемента списка инициализации
// массива функций-членов класса матрицы жесткости, предназначенных для ее вывода
#define STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT(fileFormatName, _) \
	&Matrix::WriteTo##fileFormatName,

// определелить макрос, отвечающий за развертывание формата файла матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT

// массив функций-членов класса матрицы жесткости, предназначенных для ее вывода
template
	<
		class Matrix // тип матрицы жесткости
	>
// static
const typename StiffnessMatrixIO<Matrix>::WriteMember StiffnessMatrixIO<Matrix>::_writeMembers[] =
	{
#		include "StiffnessMatrixFileFormatList.h"
	};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT


// формат файла матрицы жесткости в качестве элемента списка инициализации
// массива строк, содержащего расширения файлов матриц жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT(_, extension) \
	extension,

// определелить макрос, отвечающий за развертывание формата файла матрицы жесткости
#define STIFFNESS_MATRIX_FILE_FORMAT STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT

// массив строк, содержащий расширения файлов матриц жесткости
template
	<
		class Matrix // тип матрицы жесткости
	>
// static
const std::string StiffnessMatrixIO<Matrix>::_fileFormatsExtensions[] =
	{
#		include "StiffnessMatrixFileFormatList.h"
	};

#undef STIFFNESS_MATRIX_FILE_FORMAT

#undef STIFFNESS_MATRIX_FILE_FORMAT_INITIALIZER_ELEMENT


#endif // STIFFNESS_MATRIX_IO_HPP