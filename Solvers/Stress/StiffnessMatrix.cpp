#include "StiffnessMatrix.h"

#include "../../Fcore/Exceptions/fcExceptions.h"

#include <fstream>
#include <iomanip>


using std::ofstream;
using std::string;


// class StiffnessMatrix

// �������, ����������� � ����� �����, � ������� ������������ ��� �� �������� �������� ������� ���������
// static
const string StiffnessMatrix::_fileNamePostfix;

/**
* �����������
* @param nodesCount ���������� ����� ��������� ������������� ����
*/
StiffnessMatrix::StiffnessMatrix
	(
		size_t nodesCount
	)
	:
		_nodesCount(nodesCount)
{
	Allocate();
}

/**
* ����������� �����������
* @param otherStiffnessMatrix - ������� ���������, �� ������� ������������ �����������
*/
StiffnessMatrix::StiffnessMatrix
	(
		StiffnessMatrix&& otherStiffnessMatrix
	)
	:
		_nodesCount(otherStiffnessMatrix._nodesCount),
		_matrix(std::move(otherStiffnessMatrix._matrix))
{
	otherStiffnessMatrix._nodesCount = 0;
}

/**
* ��������� ������� ������� ���������
* @param rowIndex ������ ���� ������� ���������, ���������� �� ��������� 6 x 6
* @param rowDofIndex ������ ������ ���������� 6 x 6
* @param columnIndex ������ ������� ������� ���������, ��������� �� ��������� 6 x 6
* @param columnDofIndex ������ ������� ���������� 6 x 6
* @param elementValue �������� ������������ �������� ������� ���������
*/
void StiffnessMatrix::FillElement
	(
		size_t rowIndex,
		size_t rowDofIndex,
		size_t columnIndex,
		size_t columnDofIndex,
		double elementValue
	)
{
	_matrix[rowIndex][columnIndex][rowDofIndex][columnDofIndex] = elementValue;
}

/**
* �������� ������� 6 x 6 � ��������� �������� �����
* @param textFileStream ����� ������, � ������� ������������ ������� 6 x 6
* @param matrix6x6 - ������������ ������� 6 x 6
*/
void WriteMatrix6x6ToTextFile
	(
		ofstream& textFileStream,
		const Matrix6x6& matrix6x6
	)
{
	textFileStream << std::scientific << std::setprecision(6);
	for (const auto& row : matrix6x6)
	{
		for (const auto element : row)
		{
			textFileStream << std::setw(15) << element;
		}
		textFileStream << std::endl;
	}
	textFileStream << std::endl;
}

/**
* �������� ������� ��������� � ��������� ����
* @param textFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
*/
void StiffnessMatrix::WriteToTextFile
	(
		const string& textFileName
	)	const
{
	ofstream textFileStream(textFileName);

	if (!textFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(textFileName);
	}
	for (size_t rowIndex = 0; rowIndex < _nodesCount; ++rowIndex)
	{
		for (size_t columnIndex = 0; columnIndex < _nodesCount; ++columnIndex)
		{
			textFileStream << "(" << rowIndex << ", " << columnIndex << ")" << std::endl;
			WriteMatrix6x6ToTextFile(textFileStream, _matrix[rowIndex][columnIndex]);
		}
	}
}

/**
* �������� ������� ��������� � ��������� ���� ������� GNU Octave
* @param octaveTextFileName ������������ ���������� �����, � ������� ������������ ������ ������� ���������
*/
void StiffnessMatrix::WriteToOctaveTextFile
	(
		const string& octaveTextFileName
	)	const
{
	ofstream octaveTextFileStream(octaveTextFileName);

	if (!octaveTextFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(octaveTextFileName);
	}
	octaveTextFileStream << "[";
	for (size_t rowIndex = 0; rowIndex < _nodesCount; ++rowIndex)
	{
		for (size_t rowDofIndex = 0; rowDofIndex < FreedomsCount; ++rowDofIndex)
		{
			for (size_t columnIndex = 0; columnIndex < _nodesCount; ++columnIndex)
			{
				const Matrix6x6& matrix6x6 = _matrix[rowIndex][columnIndex];

				for (size_t columnDofIndex = 0; columnDofIndex < FreedomsCount; ++columnDofIndex)
				{
					octaveTextFileStream << matrix6x6[rowDofIndex][columnDofIndex];
					if ((columnIndex != _nodesCount - 1) || (columnDofIndex != FreedomsCount - 1))
					{
						octaveTextFileStream << ',';
					}
				}
			}
			if ((rowIndex != _nodesCount - 1) || (rowDofIndex == FreedomsCount - 1))
			{
				octaveTextFileStream << ";";
			}
		}
	}
	octaveTextFileStream << "]";
}

/**
* �������� ������� ��������� � �������� ����
* @param binaryFileName ������������ ��������� �����, � ������� ������������ ������ ������� ���������
*/
void StiffnessMatrix::WriteToBinaryFile
	(
		const string& binaryFileName
	)	const
{
	ofstream binaryFileStream(binaryFileName, std::ios::binary);

	if (!binaryFileStream.is_open())
	{
		exceptions::ThrowFileNotSaved(binaryFileName);
	}

	int matrixSize = _nodesCount * FreedomsCount;
	int nonzeroElementsCount{};

	binaryFileStream.write(reinterpret_cast<char*>(&matrixSize), sizeof(int));
	binaryFileStream.write(reinterpret_cast<char*>(&nonzeroElementsCount), sizeof(int));
	for (size_t rowIndex = 0; rowIndex < _nodesCount; ++rowIndex)
	{
		for (size_t rowDofIndex = 0; rowDofIndex < FreedomsCount; ++rowDofIndex)
		{
			for (size_t columnIndex = 0; columnIndex < _nodesCount; ++columnIndex)
			{
				const Matrix6x6& matrix6x6 = _matrix[rowIndex][columnIndex];

				for (size_t columnDofIndex = 0; columnDofIndex < FreedomsCount; ++columnDofIndex)
// TODO: �������, ����� �� ��������� ����������
#if 0
	for (size_t columnIndex = 0; columnIndex < _nodesCount; ++columnIndex)
	{
		for (size_t columnDofIndex = 0; columnDofIndex < FreedomsCount; ++columnDofIndex)
		{
			for (size_t rowIndex = 0; rowIndex < _nodesCount; ++rowIndex)
			{
				const Matrix6x6& matrix6x6 = _matrix[rowIndex][columnIndex];

				for (size_t rowDofIndex = 0; rowDofIndex < FreedomsCount; ++rowDofIndex)
#endif
				{
					const double element = matrix6x6[rowDofIndex][columnDofIndex];

					if (element != 0.0)
					{
						int fullRowIndex = FreedomsCount * rowIndex + rowDofIndex;
						int fullColumnIndex = FreedomsCount * columnIndex + columnDofIndex;

						binaryFileStream.write(reinterpret_cast<const char*>(&fullRowIndex), sizeof(int));
						binaryFileStream.write(reinterpret_cast<const char*>(&fullColumnIndex), sizeof(int));
						binaryFileStream.write(reinterpret_cast<const char*>(&element), sizeof(double));
						++nonzeroElementsCount;
					}
				}
			}
		}
	}
	binaryFileStream.seekp(sizeof(int), std::ios::beg);
	binaryFileStream.write(reinterpret_cast<const char*>(&nonzeroElementsCount), sizeof(int));
}

/**
* �������� ��������� �� ��� ������� ��������� � �������� rowIndex (��� ������� �� ��������� 6 x 6)
* @param rowIndex ������ ���� ������� ���������
* @return ��������� �� ��� ������� ��������� � �������� rowIndex
*/
PMatrix6x6& StiffnessMatrix::operator[]
	(
		size_t rowIndex
	)
{
	return _matrix[rowIndex];
}

/**
* �������� ������ ��� �������� ������� ���������
* (�� ���������� ����� ����� � �������� ������������� _nodesCount)
*/
void StiffnessMatrix::Allocate()
{
	_matrix = std::make_unique<PMatrix6x6[]>(_nodesCount);
	for (size_t rowIndex = 0; rowIndex < _nodesCount; ++rowIndex)
	{
		_matrix[rowIndex] = std::make_unique<Matrix6x6[]>(_nodesCount);
	}
}


// ����� ������������� StiffnessMatrixIO

// ��� ���� StiffnessMatrix
template class StiffnessMatrixIO<StiffnessMatrix>;