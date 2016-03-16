#include "GridSelectedArea.h"


// class GridSelectedArea

/**
* �������� ������� ����� � ���� �������������� � �������� ������� � ���������
* @param meshDataProvider - ��������� ������� � ������ �����
* @param centerPosition - ��������� ������ ���������� �������
* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
*/
void GridSelectedArea::Form
	(
		const MeshDataProvider* meshDataProvider,
		const int centerPosition[3],
		const int areaSizes[3]
	)
{
	SelectArea
		(
			meshDataProvider,
			centerPosition,
			areaSizes,
			_areaPoints,
			_subAreasPoints
		);
}

/**
* ������� ���� � ���������� �������� ����� (� ���� ���������������)
* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
* @param meshDataProvider - ��������� ������� � ������ �����
* @param centerPosition - ��������� ������ ���������� ������� � �����
* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
* @return ������� ��������� (true) ��� ����������� (false) �������� �����
*/
// static
bool GridSelectedArea::CreateFile
	(
		const std::string& selectedAreaFileName,
		const MeshDataProvider* meshDataProvider,
		const int centerPosition[3],
		const int areaSizes[3]
	)
{
	GridSelectedArea selectedArea;

	selectedArea.Form
		(
			meshDataProvider,
			centerPosition,
			areaSizes
		);

	return selectedArea.ToFile(selectedAreaFileName);
}

/**
* �������� ������ �����, �������� � ���������� �������
* @return ������ �����, �������� � ���������� �������
*/
const GridSelectedArea::PointList& GridSelectedArea::GetAreaPoints() const
{
	return _areaPoints;
}

/**
* �������� ���� � ���������� �������� ����� (� ���� ���������������)
* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
* @return ������� �������� (true) ��� ���������� (false) ������ �����
*/
bool GridSelectedArea::ToFile
	(
		const std::string& selectedAreaFileName
	)	const
{
	return WriteFile
		(
			selectedAreaFileName,
			_areaPoints,
			_subAreasPoints
		);
}

/**
* ������� ���� � ���������� �������� ����� (� ���� ���������������)
* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
* @return ������� �������� (true) ��� ���������� (false) ������ �����
*/
bool GridSelectedArea::FromFile
	(
		const std::string& selectedAreaFileName
	)
{
	return ReadFile
		(
			selectedAreaFileName,
			_areaPoints,
			_subAreasPoints
		);
}

/**
* �������� ������� ����� � ���� �������������� � �������� ������� � ���������
* @param meshDataProvider - ��������� ������� � ������ �����
* @param centerPosition - ��������� ������ ���������� �������
* @param areaSizes - ������� ���������� ������� (������ ���� ����� ���������������) � ����� �����
* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
*/
// static
void GridSelectedArea::SelectArea
	(
		const MeshDataProvider* meshDataProvider,
		const int centerPosition[3],
		const int areaSizes[3],
		PointList& selectedArea,
		PointLists& selectedSubAreas
	)
{
	selectedArea.clear();
	for (PointList& selectedSubArea : selectedSubAreas)
	{
		selectedSubArea.clear();
	}

	int areaHalfSizes[3];
	int gridSizes[3];
	int startPosition[3];
	int endPosition[3];

	meshDataProvider->GetSizes(gridSizes);
	for (int direction = 0; direction < 3; ++direction)
	{
		areaHalfSizes[direction] = areaSizes[direction] / 2;
		startPosition[direction] = std::max(0, centerPosition[direction] - areaHalfSizes[direction]);
		endPosition[direction] =
			std::min(gridSizes[direction] - 1, centerPosition[direction] + areaHalfSizes[direction]);
	}

	int position[3];

	std::copy(startPosition, startPosition + 3, position);

	bool stop = false;

	while (!stop)
	{
		int pointNumber;

		if (meshDataProvider->GetNodeByPosition(position, pointNumber))
		{
			int pointIndex = pointNumber - 1;

			selectedArea.push_back(pointIndex);
			for (int direction = 0; direction < 3; ++direction)
			{
				if (position[direction] == centerPosition[direction])
				{
					selectedSubAreas[direction].push_back(pointIndex);
				}
			}
		}

		int direction;

		for (direction = 0; direction < 3; ++direction)
		{
			position[direction]++;
			if (position[direction] < endPosition[direction])
			{
				break;
			}
			position[direction] = startPosition[direction];
		}
		if (direction == 3)
		{
			stop = true;
		}
	}
}

/**
* �������� � ���� ������ ������� ����� �����
* @return ������� �������� (true) ��� ���������� (false) ������ �����
*/
// static
bool GridSelectedArea::WritePointList
	(
		ofstream& outputFile,
		const PointList& pointList
	)
{
	if (pointList.empty())
	{
		return false;
	}

	outputFile << pointList.size() << std::endl;
	for (int pointNumber : pointList)
	{
		outputFile << pointNumber << " ";
	}
	outputFile << std::endl;

	return true;
}

/**
* ������� � ���� ������ ������� ����� �����
* @return ������� ��������� (true) ��� ����������� (false) ������ �����
*/
// static
bool GridSelectedArea::ReadPointList
	(
		ifstream& inputFile,
		PointList& pointList
	)
{
	size_t pointListSize{};
	PointList pointListTmp;

	inputFile >> pointListSize;
	if (inputFile.eof())
	{
		return false;
	}
	for (size_t index = 0; index < pointListSize; ++index)
	{
		int pointIndex{};

		inputFile >> pointIndex;
		if (inputFile.eof())
		{
			return false;
		}
		pointListTmp.push_back(pointIndex);
	}
	pointList = pointListTmp;

	return true;
}

/**
* �������� ���� � ���������� �������� ����� (� ���� ���������������)
* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
* @return ������� �������� (true) ��� ���������� (false) ������ �����
*/
// static
bool GridSelectedArea::WriteFile
	(
		const string& selectedAreaFileName,
		const PointList& selectedArea,
		const PointLists& selectedSubAreas
	)
{
	ofstream selectedAreaFile(selectedAreaFileName);

	if (!selectedAreaFile.is_open())
	{
		return false;
	}
	if (!WritePointList(selectedAreaFile, selectedArea))
	{
		return false;
	}
	for (const PointList& selectedSubArea : selectedSubAreas)
	{
		if (!WritePointList(selectedAreaFile, selectedSubArea))
		{
			return false;
		}
	}

	return true;
}

/**
* ������� ���� � ���������� �������� ����� (� ���� ���������������)
* @param selectedAreaFileName - ������������ ����� � ���������� �������� �����
* @param selectedArea - ������ ������� ����� �����, ������������ ���������� �������
* @param selectedSubAreas - ������ ������� ������� �����, ���������� ��������� ���������� �������,
* ����������������� ������������ ���� � ����������� ����� ����� ���������� �������
* @return ������� ��������� (true) ��� ����������� (false) ������ �����
*/
// static
bool GridSelectedArea::ReadFile
	(
		const std::string& selectedAreaFileName,
		PointList& selectedArea,
		PointLists& selectedSubAreas
	)
{
	ifstream selectedAreaFile(selectedAreaFileName);

	if (!selectedAreaFile.is_open())
	{
		return false;
	}

	PointList selectedAreaTmp;

	if (!ReadPointList(selectedAreaFile, selectedAreaTmp))
	{
		return false;
	}

	PointLists selectedSubAreasTmp;

	for (PointList& selectedSubAreaTmp : selectedSubAreasTmp)
	{
		if (!ReadPointList(selectedAreaFile, selectedSubAreaTmp))
		{
			return false;
		}
	}
	selectedArea = selectedAreaTmp;
	selectedSubAreas = selectedSubAreasTmp;

	return true;
}