#include "OccRectilinearGrid.h"


/**
* Функция округления (как выяснилось в VS2012 функция std::round еще не реализована)
* @param value - округляемое значение
* @return целое число, равное округленному значению
*/
template
	<
		typename Real // тип вещественного числа
	>
int round
	(
		Real value
	)
{
	return static_cast<int>(value + static_cast<Real>(0.5L));
}


OccRectilinearGrid::OccRectilinearGrid
	(
		int nl,
		int nr,
		int nel
	)
{
	nx = nl;
	ny = nr;
	nz = nel;
	for (int i = 0; i < nl; i++)
	{
		Grid.push_back(new MeshLayer(nr, nel));
	}
}

OccRectilinearGrid::OccRectilinearGrid
	(
		OccRectilinearGrid* basicGrid
	)
{
	basicGrid->GetGridSize(&nx, &ny, &nz);
	for (int i = 0; i < nx; i++)
	{
		Grid.push_back(new MeshLayer(ny, nz));
	}
	basicGrid->GetStartPoint(&startPointX, &startPointY, &startPointZ);
	basicGrid->GetElementSize(&dx, &dy, &dz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++) 
			{
				(*this)[x][y][z] = (*basicGrid)[x][y][z];
			}
		}
	}
}

OccRectilinearGrid::~OccRectilinearGrid()
{
	for (MeshLayer* layer : Grid)
	{
		delete layer;
	}
}

//вернуть размерность сетки
void OccRectilinearGrid::GetGridSize
	(
		int *x,
		int *y,
		int *z
	)	const
{
	*x = nx;
	*y = ny;
	*z = nz;
}

//вернуть начальную точку сетки
void OccRectilinearGrid::GetStartPoint
	(
		double *x,
		double *y,
		double *z
	)	const
{
	*x = startPointX;
	*y = startPointY;
	*z = startPointZ;
}

//вернуть размер элемента сетки
void OccRectilinearGrid::GetElementSize
	(
		double *x,
		double *y,
		double *z
	)	const
{
	*x = dx;
	*y = dy;
	*z = dz;
}

//вернуть количество слоев в Grid
size_t OccRectilinearGrid::GetNGrid() const
{
	return Grid.size();
}

/**
* Вернуть значение элемента сетки
* @param x - положение элемента по оси Ox
* @param y - положение элемента по оси Oy
* @param z - положение элемента по оси Oz
* @return значение элемента по данным координатам
*/
int OccRectilinearGrid::GetElement
	(
		int x,
		int y,
		int z
	)	const
{
	if (x > -1 && y > -1 && z > -1 && x < nx && y < ny && z < nz)
	{
		return GetElementInternal(x, y, z);
	}
	else
	{
		return 0;
	}
}

/**
* Получить размеры сетки (количество элменетов сетки, отложенных вдоль каждой из осей)
* @param size - массив из 3-х элементов, равных размерам сетки вдоль соответствующих осей
*/
void OccRectilinearGrid::GetSizes
	(
		int sizes[3]
	)	const
{
	sizes[0] = nx;
	sizes[1] = ny;
	sizes[2] = nz;
}

/**
* Получить данные о минимальной и максимальной координатах сетки вдоль оси Ox
* @param minCoordinate - минимальная координата сетки вдоль оси Ox
* @param maxCoordinate - максимальная координата сетки вдоль оси Ox
*/
void OccRectilinearGrid::GetMinMaxCoordinatesX
	(
		double& minCoordinate,
		double& maxCoordinate
	)	const
{
	minCoordinate = startPointX - dx / 2;
	maxCoordinate = startPointX + nx * dx - dx / 2;
}

/**
* Получить данные о минимальной и максимальной координатах сетки вдоль оси Oy
* @param minCoordinate - минимальная координата сетки вдоль оси Oy
* @param maxCoordinate - максимальная координата сетки вдоль оси Oy
*/
void OccRectilinearGrid::GetMinMaxCoordinatesY
	(
		double& minCoordinate,
		double& maxCoordinate
	)	const
{
	minCoordinate = startPointY - dy / 2;
	maxCoordinate = startPointY + ny * dy - dy / 2;
}

/**
* Получить данные о минимальной и максимальной координатах сетки вдоль оси Oz
* @param minCoordinate - минимальная координата сетки вдоль оси Oz
* @param maxCoordinate - максимальная координата сетки вдоль оси Oz
*/
void OccRectilinearGrid::GetMinMaxCoordinatesZ
	(
		double& minCoordinate,
		double& maxCoordinate
	)	const
{
	minCoordinate = startPointZ - dz / 2;
	maxCoordinate = startPointZ + nz * dz - dz / 2;
}

/**
* Получить данные о геометрических размерах сетки
* @param minCoordinates - массив из 3-х элементов, равных минимальным координатам сетки
* @param maxCoordinates - массив из 3-х элементов, равных максимальным координатам сетки
* @param elementSizes - массив из 3-х элементов, равных размерам элемента сетки вдоль соответствующих осей
*/
void OccRectilinearGrid::GetGeometricSizes
	(
		double minCoordinates[3],
		double maxCoordinates[3],
		double elementSizes[3]
	)	const
{
	GetMinMaxCoordinatesX(minCoordinates[0], maxCoordinates[0]);
	elementSizes[0] = dx;
	GetMinMaxCoordinatesY(minCoordinates[1], maxCoordinates[1]);
	elementSizes[1] = dy;
	GetMinMaxCoordinatesZ(minCoordinates[2], maxCoordinates[2]);
	elementSizes[2] = dz;
}

/**
* Получить положение узла в сетке по его номеру
* @param nodeNumber - номер узла
* @param position - найденная позиция узла в сетке
* @return признак, существует ли узел с данным номером
*/
bool OccRectilinearGrid::GetNodePosition
	(
		int nodeNumber,
		int position[3]
	)	const
{
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++) 
			{
				int currentNodeNumber = GetElementInternal(x, y, z);

				if (currentNodeNumber == nodeNumber)
				{
					position[0] = x;
					position[1] = y;
					position[2] = z;

					return true;
				}
			}
		}
	}

	return false;
}

/**
* Получить номер узла в сетке по его положению
* @param position - предполагаемая позиция узла в сетке
* @param nodeNumber - номер найденного узла
* @return признак, существует ли узел с данным положением
*/
bool OccRectilinearGrid::GetNodeByPosition
	(
		const int position[3],
		int& nodeNumber
	)	const
{
	nodeNumber = GetElementInternal(position[0], position[1], position[2]);
	if (nodeNumber == 0)
	{
		return false;
	}

	return true;
}

/**
* Получить координаты узла по его положению
* @param position - предполагаемая позиция узла в сетке
* @param coordinates - координаты найденного узла
* @return признак, существует ли узел с данным положением
*/
void OccRectilinearGrid::GetNodeCoordinatesByPosition
	(
		const int position[3],
		double coordinates[3]
	)	const
{
	coordinates[0] = startPointX + position[0] * dx;
	coordinates[1] = startPointY + position[1] * dy;
	coordinates[2] = startPointZ + position[2] * dz;
}

/**
* Получить положение узла сетки по его координатам
* @param supposedCoordinates - предполагаемые координаты узла
* @param position - положение узла в сетке
*/
void OccRectilinearGrid::GetNodePositionByCoordinates
	(
		const double supposedCoordinates[3],
		int position[3]
	)	const
{
	position[0] = static_cast<int>(round((supposedCoordinates[0] - startPointX) / dx));
	position[1] = static_cast<int>(round((supposedCoordinates[1] - startPointY) / dy));
	position[2] = static_cast<int>(round((supposedCoordinates[2] - startPointZ) / dz));
	if (position[0] >= nx)
	{
		position[0] = nx - 1;
	}
	else
	{
		if (position[0] < 0)
		{
			position[0] = 0;
		}
	}
	if (position[1] >= ny)
	{
		position[1] = ny - 1;
	}
	else
	{
		if (position[1] < 0)
		{
			position[1] = 0;
		}
	}
	if (position[2] >= nz)
	{
		position[2] = nz - 1;
	}
	else
	{
		if (position[2] < 0)
		{
			position[2] = 0;
		}
	}
}

/**
* Определить, может ли принадлежать точка с данной координатой
* ограничивающему параллелепипеду сетки
* @param direction - индекс направления оси, по которой берется координата
* @param coordinate - координата вдоль оси
*/
bool OccRectilinearGrid::IsCorrectCoordinate
	(
		int direction,
		double coordinate
	)	const
{
	double minCoordinate = 0.0;
	double maxCoordinate = 0.0;

	switch (direction)
	{
	case 0:
		GetMinMaxCoordinatesX(minCoordinate, maxCoordinate);
		break;
		
	case 1:
		GetMinMaxCoordinatesY(minCoordinate, maxCoordinate);
		break;
		
	case 2:
		GetMinMaxCoordinatesZ(minCoordinate, maxCoordinate);
		break;
	}

	return (coordinate >= minCoordinate) && (coordinate <= maxCoordinate);
}

// Перегрузка оператора []
const MeshLayer& OccRectilinearGrid::operator[]
	(
		int index
	)	const
{
	return *Grid[index];
}

// Перегрузка оператора []
MeshLayer& OccRectilinearGrid::operator[]
	(
		int index
	)
{
	return *Grid[index];
}

//установить начальную точку сетки
void OccRectilinearGrid::SetStartPoint
	(
		double x,
		double y,
		double z
	)
{
	startPointX = x;
	startPointY = y;
	startPointZ = z;
}

//установить размер элемента сетки
void OccRectilinearGrid::SetElementSize
	(
		double x,
		double y,
		double z
	)
{
	dx = x;
	dy = y;
	dz = z;
}

//Добавить слой
void OccRectilinearGrid::AddLayer
	(
		MeshLayer *layer
	)
{
	Grid.push_back(layer);
}

/**
* Вернуть значение элемента сетки
* (без проверок на попадание в ограничивающий параллелепипед)
* @param x - положение элемента по оси Ox
* @param y - положение элемента по оси Oy
* @param z - положение элемента по оси Oz
* @return значение элемента по данным координатам
*/
int OccRectilinearGrid::GetElementInternal
	(
		int x,
		int y,
		int z
	)	const
{
	return (*Grid.at(x))[y][z];
}

int OccRectilinearGrid::EnumerateMesh()
{
	int index = 1;
	int nx, ny, nz, nodesCount;

	GetGridSize(&nx, &ny, &nz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++) 
			{
				int tmp = (*this)[x][y][z];

				if (tmp > 0)
				{
					(*this)[x][y][z] = index;
					index++;
				}
			}
		}
	}
	nodesCount = index - 1;

	return nodesCount;
}

int OccRectilinearGrid::DeenumerateMesh()
{
	int index = 1;
	int nx, ny, nz, nodesCount;

	GetGridSize(&nx, &ny, &nz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++) 
			{
				int tmp = (*this)[x][y][z];

				if (tmp > 0)
				{
					(*this)[x][y][z] = 1;
					index++;
				}
			}
		}
	}
	nodesCount = index - 1;

	return nodesCount;
}