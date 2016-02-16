#include "StiffnessRHS.h"

#include "../../Fcore/Exceptions/fcExceptions.h"

#include <fstream>
#include <iomanip>


using std::ifstream;
using std::ofstream;


/**
* ��������� ������ ������� ������ ������������ �������� ������� ��� ��������� ��������� ������������� ����
* � ������� ������ ������ ���� �� ������ ������� ��������� (��������������� �������� ������� ����� ��������)
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
* @param isSealedFlagsList - ������ ������� ������ ������������ �������� �������
* ��������� ��������� �������������
*/
void ApplySealedElementsToStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector,
		const IsSealedFlagsList& isSealedFlagsList
	)
{
	if (stiffnessRHSVector.size() != isSealedFlagsList.size() * FreedomsCount)
	{
		exceptions::ThrowMessage
			(
				"StressStrainFortranIterativeSolver::ApplySealedElementsToStiffnessRHSVector:\n"
				"Sizes of stiffnessRHSVector and isSealedFlagsList are incompatible!"
			);
	}

	int elementIndex{};

	for (const IsSealedFlags& isSealedFlags : isSealedFlagsList)
	{
		for (int dofIndex = 0; dofIndex < 6; ++dofIndex)
		{
			if (isSealedFlags[dofIndex])
			{
				stiffnessRHSVector[elementIndex] = 0.0;
			}
			++elementIndex;
		}
	}
}

/**
* ������������� ������ ������ ������ ���� �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void InvertStiffnessRHSVector
	(
		StiffnessRHSVector& stiffnessRHSVector
	)
{
	for (double& element : stiffnessRHSVector)
	{
		if (element != 0.0)
		{
			element = -element;
		}
	}
}

/**
* �������� ������ ������ ������ ���� �� ������ ������� ��������� � ��������� ����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void WriteStiffnessRHSVectorToTextFile
	(
		const string& stiffnessRHSTextFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	)
{
	ofstream stiffnessRHSTextFile(stiffnessRHSTextFileName);

	stiffnessRHSTextFile << std::scientific << std::setprecision(6);

	int elementIndex{};

	for (double element : stiffnessRHSVector)
	{
		if (elementIndex % 6 == 0)
		{
			stiffnessRHSTextFile << "Node index = " << elementIndex / 6 << std::endl;
		}
		stiffnessRHSTextFile << std::setw(15) << element << std::endl;
	}
}

/**
* �������� ������ ������ ������ ���� �� ������ ������� ��������� � �������� ����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @param stiffnessRHSVector - ������ ������ ������ ���� �� ������ ������� ���������
*/
void WriteStiffnessRHSVectorToBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName,
		const StiffnessRHSVector& stiffnessRHSVector
	)
{
	ofstream stiffnessRHSBinaryFile(stiffnessRHSBinaryFileName, std::ios::binary);

	if (stiffnessRHSBinaryFile.is_open())
	{
		stiffnessRHSBinaryFile.write
			(
				reinterpret_cast<const char*>(stiffnessRHSVector.data()),
				stiffnessRHSVector.size() * sizeof(double)
			);
	}
}

/**
* ������� ������ ������ ������ ���� �� ������ ������� ��������� �� ��������� �����
* @param stiffnessRHSTextFileName - ������������ ����� ��� ������ ������� ������ ������ ����
* �� ������ ������� ���������
* @return ��������� ������ ������ ������ ���� �� ������ ������� ���������
*/
StiffnessRHSVector ReadStiffnessRHSVectorFromBinaryFile
	(
		const std::string& stiffnessRHSBinaryFileName
	)
{
	ifstream stiffnessRHSBinaryFile(stiffnessRHSBinaryFileName, std::ios::binary | std::ios::ate);

	if (!stiffnessRHSBinaryFile.is_open())
	{
		// TODO: �������� ������ ����������
		return {};
	}

	auto stiffnessRHSBinaryFileSize = stiffnessRHSBinaryFile.tellg();

	if (stiffnessRHSBinaryFileSize <= 0)
	{
		// TODO: �������� ������ ����������
		return {};
	}
	if (stiffnessRHSBinaryFileSize % sizeof(double) != 0)
	{
		// TODO: �������� ������ ����������
		return {};
	}
	stiffnessRHSBinaryFile.seekg(0, std::ios::beg);

	const size_t stiffnessRHSVectorSize = static_cast<size_t>(stiffnessRHSBinaryFileSize) / sizeof(double);
	StiffnessRHSVector stiffnessRHSVector(stiffnessRHSVectorSize);

	stiffnessRHSBinaryFile.read
		(
			reinterpret_cast<char*>(stiffnessRHSVector.data()),
			stiffnessRHSVector.size() * sizeof(double)
		);

	return stiffnessRHSVector;
}