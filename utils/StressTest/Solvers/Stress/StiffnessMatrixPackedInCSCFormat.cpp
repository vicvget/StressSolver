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

// �������, ����������� � ����� �����, � ������� ������������ ��� �� �������� �������� ������� ���������
// static
const std::string StiffnessMatrixPackedInCSCFormat::_fileNamePostfix = "PackedCSCF";

/**
* �����������
* @param nodesCount ���������� ����� ��������� ������������� ����
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
* �������� ������ ������ (�������������) ������� ���������
* @return ������ ������ (�������������) ������� ���������
*/
size_t StiffnessMatrixPackedInCSCFormat::GetMatrixSize() const
{
	return _matrixSize;
}

/**
* �������� ���������� ��������� ��������� � ������� ���������
* @return ���������� ��������� ��������� � ������� ���������
*/
size_t StiffnessMatrixPackedInCSCFormat::GetMatrixNonZeroElementsCount() const
{
	return _values.size();
}

/**
* �������� ��������� �� ������ ��������� ������� ���������
* @return ��������� �� ������ ��������� ������� ���������
*/
const double* StiffnessMatrixPackedInCSCFormat::GetElements() const
{
	return _values.data();
}

/**
* ����������� �������� ������� ��������� � ������ ������
* @param elements - ������, � ������� ������������ ����������� ��������� ������� ���������
*/
void StiffnessMatrixPackedInCSCFormat::CopyElements
	(
		double* elements
	)	const
{
	std::copy(_values.begin(), _values.end(), elements);
}

/**
* ����������� ������� ����� ��������� � ������� ��������� � ������ ������
* @param rowIds - ������, � ������� ������������ ����������� �������� ����� ��������� � ������� ���������
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
* ����������� ������� ������ ��������� ���������� ������� � ���������� ������� � ������ ������
* @param columnsFirstIds - ������, � ������� ������������ ����������� �������� ������ ���������
* ���������� ������� � ���������� ������� � ������ ������
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
* ��������� ������� ������� ���������
* @param rowIndex ������ ������ ������� ���������
* @param columnIndex ������ ������� ������� ���������
* @param elementValue �������� ������������ �������� ������� ���������
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
* ��������� ������� ������� ���������
* @param rowIndex ������ ���� ������� ���������, ���������� �� ��������� 6 x 6
* @param rowDofIndex ������ ������ ���������� 6 x 6
* @param columnIndex ������ ������� ������� ���������, ��������� �� ��������� 6 x 6
* @param columnDofIndex ������ ������� ���������� 6 x 6
* @param elementValue �������� ������������ �������� ������� ���������
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
* �������� ������� ��������� � ��������� ����
* @param textFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
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
* �������� ������� ��������� � ��������� ���� ������� GNU Octave
* @param octaveTextFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
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
* �������� ������� ��������� � �������� ����
* @param binaryFileName ������������ ��������� �����, � ������� ������������ ������ ������� ���������
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
* ������� ������� ��������� �� �������� ����
* @param binaryFileName ������������ ��������� �����, �� �������� ������������ ���������� ������� ���������
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


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrixPackedInCSCFormat
template class StiffnessMatrixIO<StiffnessMatrixPackedInCSCFormat>;