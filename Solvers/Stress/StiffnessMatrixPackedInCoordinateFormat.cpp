#include "StiffnessMatrixPackedInCoordinateFormat.h"

#include "AuxiliaryStressStuff.h"
#include "../../Fcore/Exceptions/fcExceptions.h"

#include <fstream>


using std::ofstream;
using std::ifstream;
using std::string;


// class StiffnessMatrixPackedInCoordinateFormat

// �������, ����������� � ����� �����, � ������� ������������ ��� �� �������� �������� ������� ���������
// static
const std::string StiffnessMatrixPackedInCoordinateFormat::_fileNamePostfix = "PackedCF";

/**
* �����������
* @param nodesCount ���������� ����� ��������� ������������� ����
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
* ��������� ������� ������� ���������
* @param rowIndex ������ ������ ������� ���������
* @param columnIndex ������ ������� ������� ���������
* @param elementValue �������� ������������ �������� ������� ���������
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
* ��������� ������� ������� ���������
* @param rowIndex ������ ���� ������� ���������, ���������� �� ��������� 6 x 6
* @param rowDofIndex ������ ������ ���������� 6 x 6
* @param columnIndex ������ ������� ������� ���������, ��������� �� ��������� 6 x 6
* @param columnDofIndex ������ ������� ���������� 6 x 6
* @param elementValue �������� ������������ �������� ������� ���������
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
* �������� ������� ��������� � ��������� ����
* @param textFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
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
* �������� ������� ��������� � ��������� ���� ������� GNU Octave
* @param octaveTextFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
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
* �������� ������� ��������� � �������� ����
* @param binaryFileName ������������ ��������� �����, � ������� ������������ ������ ������� ���������
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
* ������� ������� ��������� �� �������� ����
* @param binaryFileName ������������ ��������� �����, �� �������� ������������ ���������� ������� ���������
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


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrixPackedInCoordinateFormat
template class StiffnessMatrixIO<StiffnessMatrixPackedInCoordinateFormat>;