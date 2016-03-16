#include "GridSelectedArea.h"


// class GridSelectedArea

/**
* Выделить область сетки в виде паралелепипеда с заданным центром и размерами
* @param meshDataProvider - провайдер доступа к данным сетки
* @param centerPosition - положение центра выделяемой области
* @param areaSizes - размеры выделяемой области (массив длин ребер параллелепипеда) в узлах сетки
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
* Создать файл с выделенной областью сетки (в виде параллелепипеда)
* @param selectedAreaFileName - наименование файла с выделенной областью сетки
* @param meshDataProvider - провайдер доступа к данным сетки
* @param centerPosition - положение центра выделенной области в сетке
* @param areaSizes - размеры выделенной области (массив длин ребер параллелепипеда) в узлах сетки
* @return признак успешного (true) или неуспешного (false) создания файла
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
* Получить список узлов, входящих в выделенную область
* @return список узлов, входящих в выделенную область
*/
const GridSelectedArea::PointList& GridSelectedArea::GetAreaPoints() const
{
	return _areaPoints;
}

/**
* Записать файл с выделенной областью сетки (в виде параллелепипеда)
* @param selectedAreaFileName - наименование файла с выделенной областью сетки
* @return признак успешной (true) или неуспешной (false) записи файла
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
* Считать файл с выделенной областью сетки (в виде параллелепипеда)
* @param selectedAreaFileName - наименование файла с выделенной областью сетки
* @return признак успешной (true) или неуспешной (false) записи файла
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
* Выделить область сетки в виде паралелепипеда с заданным центром и размерами
* @param meshDataProvider - провайдер доступа к данным сетки
* @param centerPosition - положение центра выделяемой области
* @param areaSizes - размеры выделяемой области (массив длин ребер параллелепипеда) в узлах сетки
* @param selectedArea - список номеров точек сетки, составляющих выделяемую область
* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделяемой области,
* перпендикулярными координатным осям и проходящими через центр выделяемой области
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
* Записать в файл список номеров узлов сетки
* @return признак успешной (true) или неуспешной (false) записи файла
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
* Считать в файл список номеров узлов сетки
* @return признак успешного (true) или неуспешного (false) чтения файла
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
* Записать файл с выделенной областью сетки (в виде параллелепипеда)
* @param selectedAreaFileName - наименование файла с выделенной областью сетки
* @param selectedArea - список номеров точек сетки, составляющих выделенную область
* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделенной области,
* перпендикулярными координатным осям и проходящими через центр выделенной области
* @return признак успешной (true) или неуспешной (false) записи файла
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
* Считать файл с выделенной областью сетки (в виде параллелепипеда)
* @param selectedAreaFileName - наименование файла с выделенной областью сетки
* @param selectedArea - список номеров точек сетки, составляющих выделенную область
* @param selectedSubAreas - массив списков номеров точек, являющихся сечениями выделенной области,
* перпендикулярными координатным осям и проходящими через центр выделенной области
* @return признак успешного (true) или неуспешного (false) чтения файла
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