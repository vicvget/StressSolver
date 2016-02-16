#include "StiffnessMatrixPackedInCoordinateFormat.h"

#include "AuxiliaryStressStuff.h"
#include "../../Fcore/Exceptions/fcExceptions.h"

#include <fstream>


using std::ofstream;
using std::ifstream;
using std::string;


// class StiffnessMatrixPackedInCoordinateFormat

// суффикс, добавляемый к имени файла, в который записывается или из которого читается матрица жесткости
// static
const std::string StiffnessMatrixPackedInCoordinateFormat::_fileNamePostfix = "PackedCF";

/**
* Конструктор
* @param nodesCount количество узлов сеточного представления тела
*/
StiffnessMatrixPackedInCoordinateFormat::StiffnessMatrixPackedInCoordinateFormat
	(
		size_t nodesCount
	)
	:
		_matrixSize(nodesCount * FreedomsCount)
{

}

/**
* Заполнить элемент матрицы жесткости
* @param rowIndex индекс строки матрицы жесткости
* @param columnIndex индекс столбца матрицы жесткости
* @param elementValue значение заполняемого элемента матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::FillElement
	(
		size_t rowIndex,
		size_t columnIndex,
		double elementValue
	)
{
	_matrixElements.push_back({static_cast<int>(rowIndex), static_cast<int>(columnIndex), elementValue});
}

/**
* Заполнить элемент матрицы жесткости
* @param rowIndex индекс ряда матрицы жесткости, состоящего из подматриц 6 x 6
* @param rowDofIndex индекс строки подматрицы 6 x 6
* @param columnIndex индекс колонки матрицы жесткости, состоящей из подматриц 6 x 6
* @param columnDofIndex индекс столбца подматрицы 6 x 6
* @param elementValue значение заполняемого элемента матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::FillElement
	(
		size_t rowIndex,
		size_t rowDofIndex,
		size_t columnIndex,
		size_t columnDofIndex,
		double elementValue
	)
{
	FillElement
		(
			rowIndex * FreedomsCount + rowDofIndex,
			columnIndex * FreedomsCount + columnDofIndex,
			elementValue
		);
}

/**
* Записать матрицу жесткости в текстовый файл
* @param textFileName наименование текстового файла, в который производится запись матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::WriteToTextFile
	(
		const std::string& textFileName
	)	const
{
	ofstream textFileStream(textFileName);

	if (!textFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(textFileName);
	}
	textFileStream << "Matrix size: " << _matrixSize << std::endl;
	textFileStream << "Nonzero elements count: " << _matrixElements.size() << std::endl;
	for (const MatrixElement& matrixElement : _matrixElements)
	{
		textFileStream
			<< '('
				<< matrixElement.rowIndex << ','
				<< matrixElement.columnIndex
			<< "): "
			<< matrixElement.elementValue
			<< std::endl;
	}
}

/**
* Записать матрицу жесткости в текстовый файл формата GNU Octave
* @param octaveTextFileName наименование текстового файла, в который производится запись матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::WriteToOctaveTextFile
	(
		const std::string& octaveTextFileName
	)	const
{
	exceptions::ThrowMessage
		(
			"Write to Octave text file is not implemented in class StiffnessMatrixPackedInCoordinateFormat!"
		);
}

/**
* Записать матрицу жесткости в двоичный файл
* @param binaryFileName наименование двоичного файла, в который производится запись матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::WriteToBinaryFile
	(
		const std::string& binaryFileName
	)	const
{
	ofstream binaryFileStream(binaryFileName, ofstream::binary);

	if (!binaryFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(binaryFileName);
	}

	const int matrixSize = static_cast<int>(_matrixSize);
	const int nonzeroElementsCount = static_cast<int>(_matrixElements.size());

	binaryFileStream.write(reinterpret_cast<const char*>(&matrixSize), sizeof(int));
	binaryFileStream.write(reinterpret_cast<const char*>(&nonzeroElementsCount), sizeof(int));
	if (!_matrixElements.empty())
	{
		binaryFileStream.write
			(
				reinterpret_cast<const char*>(_matrixElements.data()),
				_matrixElements.size() * sizeof(MatrixElement)
			);
	}
}

/**
* Считать матрицу жесткости из двоичный файл
* @param binaryFileName наименование двоичного файла, из которого производится считывание матрицы жесткости
*/
void StiffnessMatrixPackedInCoordinateFormat::ReadFromBinaryFile
	(
		const std::string& binaryFileName
	)
{
	ifstream binaryFileStream(binaryFileName, ifstream::binary);

	if (!binaryFileStream.is_open())
	{
		exceptions::ThrowFileNotOpened(binaryFileName);
	}

	int matrixSize{};
	int nonzeroElementsCount{};

	binaryFileStream.read(reinterpret_cast<char*>(&matrixSize), sizeof(int));
	_matrixSize = matrixSize;
	binaryFileStream.read(reinterpret_cast<char*>(&nonzeroElementsCount), sizeof(int));
	_matrixElements.resize(nonzeroElementsCount);
	binaryFileStream.read
		(
			reinterpret_cast<char*>(_matrixElements.data()),
			_matrixElements.size() * sizeof(MatrixElement)
		);
}


// Явная конкретизация StiffnessMatrixIO

// для типа StiffnessMatrixPackedInCoordinateFormat
template class StiffnessMatrixIO<StiffnessMatrixPackedInCoordinateFormat>;