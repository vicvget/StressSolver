#include "StiffnessMatrixPackedInCSCFormat.h"

#include "AuxiliaryStressStuff.h"
#include "StiffnessMatrixPackedInCoordinateFormat.h"
#include "../../Fcore/Exceptions/fcExceptions.h"

#include <algorithm>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::string;


// class StiffnessMatrixPackedInCSCFormat

// суффикс, добавляемый к имени файла, в который записывается или из которого читается матрица жесткости
// static
const std::string StiffnessMatrixPackedInCSCFormat::_fileNamePostfix = "PackedCSCF";

/**
* Конструктор
* @param nodesCount количество узлов сеточного представления тела
*/
StiffnessMatrixPackedInCSCFormat::StiffnessMatrixPackedInCSCFormat
	(
		size_t nodesCount
	)
	:
		_matrixSize(nodesCount * FreedomsCount)
{
	_columnBeginnings.reserve(_matrixSize + 1);
}

/**
* Получить размер полной (неупакованной) матрицы жесткости
* @return размер полной (неупакованной) матрицы жесткости
*/
size_t StiffnessMatrixPackedInCSCFormat::GetMatrixSize() const
{
	return _matrixSize;
}

/**
* Получить количество ненулевых элементов в матрице жесткости
* @return количество ненулевых элементов в матрице жесткости
*/
size_t StiffnessMatrixPackedInCSCFormat::GetMatrixNonZeroElementsCount() const
{
	return _values.size();
}

/**
* Получить указатель на массив элементов матрицы жесткости
* @return указатель на массив элементов матрицы жесткости
*/
const double* StiffnessMatrixPackedInCSCFormat::GetElements() const
{
	return _values.data();
}

/**
* Скопировать элементы матрицы жесткости в другой массив
* @param elements - массив, в который производится копирвоание элементов матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::CopyElements
	(
		double* elements
	)	const
{
	std::copy(_values.begin(), _values.end(), elements);
}

/**
* Скопировать индексы строк элементов в матрице жесткости в другой массив
* @param rowIds - массив, в который производится копирование индексов строк элементов в матрице жесткости
*/
void StiffnessMatrixPackedInCSCFormat::CopyRowIds
	(
		int* rowIds
	)	const
{
	std::transform
		(
			_rows.begin(),
			_rows.end(),
			rowIds,
			[](size_t element)
			{
				return static_cast<int>(element);
			}
		);
}

/**
* Скопировать индексы первых элементов следующего столбца в предыдущем массиве в другой массив
* @param columnsFirstIds - массив, в который производится копирвоание индексов первых элементов
* следующего столбца в предыдущем массиве в другой массив
*/
void StiffnessMatrixPackedInCSCFormat::CopyColumnFirstElements
	(
		int* columnsFirstIds
	)	const
{
	std::transform
		(
			_columnBeginnings.begin(),
			_columnBeginnings.end(),
			columnsFirstIds,
			[](size_t element)
			{
				return static_cast<int>(element);
			}
		);
}

/**
* Заполнить элемент матрицы жесткости
* @param rowIndex индекс строки матрицы жесткости
* @param columnIndex индекс столбца матрицы жесткости
* @param elementValue значение заполняемого элемента матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::FillElement
	(
		size_t rowIndex,
		size_t columnIndex,
		double elementValue
	)
{
	if (columnIndex + 2 <= _columnBeginnings.size())
	{
		const size_t currentColumnBeginning = _columnBeginnings[columnIndex];
		const size_t nextColumnBeginning = _columnBeginnings[columnIndex + 1];
		const auto columnBegin = _rows.begin() + currentColumnBeginning;
		const auto columnEnd = _rows.begin() + nextColumnBeginning;
		const auto lowerBound = std::lower_bound(columnBegin, columnEnd, rowIndex);
		const size_t lowerBoundIndex = lowerBound - _rows.begin();

		if ((lowerBound != columnEnd) && (*lowerBound == rowIndex))
		{
			_values[lowerBoundIndex] = elementValue;
		}
		else
		{
			_rows.insert(lowerBound, 1, rowIndex);
			_values.insert(_values.begin() + lowerBoundIndex, 1, elementValue);
			std::for_each
				(
					_columnBeginnings.begin() + columnIndex + 1,
					_columnBeginnings.end(),
					[](size_t& columnBeginning)
					{
						++columnBeginning;
					}
				);
		}
	}
	else
	{
		_columnBeginnings.insert
			(
				_columnBeginnings.end(),
				columnIndex + 1 - _columnBeginnings.size(),
				_rows.size()
			);
		_columnBeginnings.push_back(_rows.size() + 1);
		_rows.push_back(rowIndex);
		_values.push_back(elementValue);
	}
}

/**
* Заполнить элемент матрицы жесткости
* @param rowIndex индекс ряда матрицы жесткости, состоящего из подматриц 6 x 6
* @param rowDofIndex индекс строки подматрицы 6 x 6
* @param columnIndex индекс колонки матрицы жесткости, состоящей из подматриц 6 x 6
* @param columnDofIndex индекс столбца подматрицы 6 x 6
* @param elementValue значение заполняемого элемента матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::FillElement
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
void StiffnessMatrixPackedInCSCFormat::WriteToTextFile
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
	textFileStream << "Nonzero elements count: " << _values.size() << std::endl;
	for
		(
			size_t columnBeginningIndex = 0;
			columnBeginningIndex < _columnBeginnings.size() - 1;
			++columnBeginningIndex
		)
	{
		const size_t currentColumnBeginning = _columnBeginnings[columnBeginningIndex];
		const size_t nextColumnBeginning = _columnBeginnings[columnBeginningIndex + 1];
		const size_t elementsInColumnCount = nextColumnBeginning - currentColumnBeginning;

		for (size_t rowIndex = currentColumnBeginning; rowIndex < nextColumnBeginning; ++rowIndex)
		{
			textFileStream
				<< '('
					<< _rows[rowIndex] << ','
					<< columnBeginningIndex
				<< "): "
				<< _values[rowIndex]
				<< std::endl;
		}
	}
}

/**
* Записать матрицу жесткости в текстовый файл формата GNU Octave
* @param octaveTextFileName наименование текстового файла, в который производится запись матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::WriteToOctaveTextFile
	(
		const std::string& octaveTextFileName
	)	const
{
	exceptions::ThrowMessage
		(
			"Write to Octave text file is not implemented in class StiffnessMatrixPackedInCSCFormat!"
		);
}

/**
* Записать матрицу жесткости в двоичный файл
* @param binaryFileName наименование двоичного файла, в который производится запись матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::WriteToBinaryFile
	(
		const std::string& binaryFileName
	)	const
{
	ofstream binaryFileStream(binaryFileName, ofstream::binary);

	if (!binaryFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(binaryFileName);
	}

	const int matrixSize = _matrixSize;
	const int nonzeroElementsCount = _values.size();

	binaryFileStream.write(reinterpret_cast<const char*>(&matrixSize), sizeof(int));
	binaryFileStream.write(reinterpret_cast<const char*>(&nonzeroElementsCount), sizeof(int));
	for
		(
			size_t columnBeginningIndex = 0;
			columnBeginningIndex < _columnBeginnings.size() - 1;
			++columnBeginningIndex
		)
	{
		const size_t currentColumnBeginning = _columnBeginnings[columnBeginningIndex];
		const size_t nextColumnBeginning = _columnBeginnings[columnBeginningIndex + 1];
		const size_t elementsInColumnCount = nextColumnBeginning - currentColumnBeginning;

		for (size_t rowIndex = currentColumnBeginning; rowIndex < nextColumnBeginning; ++rowIndex)
		{
			MatrixElement matrixElement{_rows[rowIndex], columnBeginningIndex, _values[rowIndex]};

			binaryFileStream.write(reinterpret_cast<const char*>(&matrixElement), sizeof(MatrixElement));
		}
	}
}

/**
* Считать матрицу жесткости из двоичный файл
* @param binaryFileName наименование двоичного файла, из которого производится считывание матрицы жесткости
*/
void StiffnessMatrixPackedInCSCFormat::ReadFromBinaryFile
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
	_values.reserve(nonzeroElementsCount);
	_rows.reserve(nonzeroElementsCount);

	for (int nonZeroElementIndex = 0; nonZeroElementIndex < nonzeroElementsCount; ++nonZeroElementIndex)
	{
		MatrixElement matrixElement;

		binaryFileStream.read
			(
				reinterpret_cast<char*>(&matrixElement),
				sizeof(MatrixElement)
			);
		FillElement(matrixElement.rowIndex, matrixElement.columnIndex, matrixElement.elementValue);
	}
}


// Явная конкретизация StiffnessMatrixIO

// для типа StiffnessMatrixPackedInCSCFormat
template class StiffnessMatrixIO<StiffnessMatrixPackedInCSCFormat>;