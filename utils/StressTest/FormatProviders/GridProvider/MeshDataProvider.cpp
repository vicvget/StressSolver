#include "MeshDataProvider.h"
#include "RLCControlReader.h"
#include <stdexcept>

using std::ios;


// class MeshDataProvider

MeshDataProvider::MeshDataProvider
	(
		OccRectilinearGrid* newgrid
	)
	:
		_grid(newgrid)
{
	_grid->GetElementSize(&_dx, &_dy, &_dz);
}

MeshDataProvider::MeshDataProvider
	(
		OccRectilinearGrid* newgrid,
		const vector<BoundaryCondition*>& surfaces
	)
	:
		_grid(newgrid),
		_surfaces(surfaces)
{
	_grid->GetElementSize(&_dx, &_dy, &_dz);
}

MeshDataProvider::MeshDataProvider
	(
		OccRectilinearGrid* newgrid,
		const vector<BoundaryCondition*>& surfaces,
		const vector<BoundaryNormal*>& normals
	)
	:
		_grid(newgrid),
		_surfaces(surfaces),
		_normals(normals)
{
	_grid->GetElementSize(&_dx, &_dy, &_dz);
}

MeshDataProvider::MeshDataProvider
	(
		const string& gridFile
	)
{
	if (RLCHeaderFileProvider::IsNewRLCFormat(gridFile))
	{
		RLCControlReader tmpreader;
		tmpreader.ReadMeshFromFile(gridFile.c_str());
		_grid = tmpreader.GetGrid();
		_surfaces = tmpreader.GetBC();
		_normals = tmpreader.GetNormals();
		//_freeSolverGridParams = tmpreader.GetFreeSolverGridParams();
		if (_grid != nullptr)
		{
			_grid->GetElementSize(&_dx, &_dy, &_dz);
			_nodesCount = tmpreader.GetCount();
		}
	}
	else
	{
		ifstream ifs(gridFile, ios::binary);
		if (ifs.is_open())
		{
			RLCControlReader tmpreader;

			_grid = tmpreader.fillGrid(&ifs);
			//std::cout << "######################" << std::endl;
			//std::cout << "FILL GRID " << std::endl;
			_surfaces = tmpreader.GetBC();
			//std::cout << "######################" << std::endl;
			//std::cout << "GET BC " << std::endl;
			if (_grid != nullptr)
			{
				//std::cout << "######################" << std::endl;
				//std::cout << "GRID NOT NULL: "<< std::endl;
				_grid->GetElementSize(&_dx, &_dy, &_dz);
				_nodesCount = tmpreader.GetCount();
				// TODO: Debug
			}
		}
		ifs.close();
	}
}

//
OccRectilinearGrid* MeshDataProvider::GetGrid()
{
	return _grid;
}

//
const vector<BoundaryCondition*>& MeshDataProvider::GetBC() const
{
	return _surfaces;

}

const vector<BoundaryNormal*>& MeshDataProvider::GetNormals() const
{
	return _normals;
}

//FreeSolverGridParams* MeshDataProvider::GetFreeSolverGridParams() const
//{
//	return _freeSolverGridParams;
//}

//
bool MeshDataProvider::IsGridLoaded() const
{
	return _grid != nullptr;
}

//
int MeshDataProvider::GetCount() const
{
	return _nodesCount;
}

//
int MeshDataProvider::GetNodesNumber() const
{
	return _nodesCount;
}

/**
* Получить количество точек в ограничивающем сетку параллелепипеде
* @return количество точек в ограничивающем сетку параллелепипеде
*/
int MeshDataProvider::GetPointsCountInBox() const
{
	if (_grid == NULL)
	{
		return 0;
	}

	int xSize = 0;
	int ySize = 0;
	int zSize = 0;

	GetGridSize
		(
			xSize,
			ySize,
			zSize
		);

	return xSize * ySize * zSize;
}

/**
* Получить количество граничных поверхностей в сетке
* @return количество граничных поверхностей в сетке
*/
size_t MeshDataProvider::GetBondaryConditionsCount() const
{
	return _surfaces.size();
}

/**
* Получить размеры сетки в узлах по каждому из координатных направлений
* @param xSize - размер сетки вдоль оси Ox
* @param ySize - размер сетки вдоль оси Oy
* @param zSize - размер сетки вдоль оси Oz
*/
void MeshDataProvider::GetGridSize
	(
		int& xSize,
		int& ySize,
		int& zSize
	)	const
{
	if (_grid == NULL)
	{
		return;
	}
	_grid->GetGridSize
		(
			&xSize,
			&ySize,
			&zSize
		);
}

/**
* Получить координаты начальной точки сетки
* @param xStartPoint - координата начальной точки по оси Ox
* @param yStartPoint - координата начальной точки по оси Oy
* @param zStartPoint - координата начальной точки по оси Oz
*/
void MeshDataProvider::GetStartPoint
	(
		double& xStartPoint,
		double& yStartPoint,
		double& zStartPoint
	)	const
{
	if (_grid == NULL)
	{
		return;
	}
	_grid->GetStartPoint
		(
			&xStartPoint,
			&yStartPoint,
			&zStartPoint
		);
}

/**
* Получить размер элемента (узла) сетки по каждому из координатных направлений
* @param xElementSize - размер элемента (узла) сетки вдоль оси Ox
* @param yElementSize - размер элемента (узла) сетки вдоль оси Oy
* @param zElementSize - размер элемента (узла) сетки вдоль оси Oz
*/
void MeshDataProvider::GetElementSize
	(
		double& xElementSize,
		double& yElementSize,
		double& zElementSize
	)	const
{
	if (_grid == NULL)
	{
		return;
	}
	_grid->GetElementSize
		(
			&xElementSize,
			&yElementSize,
			&zElementSize
		);
}

/**
* Получить номер элемента сетки с данными координатами
* @param xCoordinate - координата элемента по оси Ox
* @param yCoordinate - координата элемента по оси Oy
* @param zCoordinate - координата элемента по оси Oz
*/
int MeshDataProvider::GetPoint
	(
		int xCoordinate,
		int yCoordinate,
		int zCoordinate
	)	const
{
	if (_grid == NULL)
	{
		return -1;
	}

	return _grid->GetElement
		(
			xCoordinate,
			yCoordinate,
			zCoordinate
		);
}

/**
* Получить размеры сетки (количество элменетов сетки, отложенных вдоль каждой из осей)
* @param size - массив из 3-х элементов, равных размерам сетки вдоль соответствующих осей
*/
void MeshDataProvider::GetSizes
	(
		int sizes[3]
	)	const
{
	if (_grid == nullptr)
	{
		return;
	}

	return _grid->GetSizes(sizes);
}

/**
* Получить данные о геометрических размерах сетки
* @param minCoordinates - массив из 3-х элементов, равных минимальным координатам сетки
* @param maxCoordinates - массив из 3-х элементов, равных максимальным координатам сетки
* @param elementSizes - массив из 3-х элементов, равных размерам элемента сетки вдоль соответствующих осей
*/
void MeshDataProvider::GetGeometricSizes
	(
		double minCoordinates[3],
		double maxCoordinates[3],
		double elementSizes[3]
	)	const
{
	if (_grid == nullptr)
	{
		return;
	}

	return _grid->GetGeometricSizes(minCoordinates, maxCoordinates, elementSizes);
}

/**
* Получить положение узла в сетке по его номеру
* @param nodeNumber - номер узла
* @param position - найденная позиция узла в сетке
* @return признак, существует ли узел с данным номером
*/
bool MeshDataProvider::GetNodePosition
	(
		int nodeNumber,
		int position[3]
	)	const
{
	if (_grid == nullptr)
	{
		return false;
	}

	return _grid->GetNodePosition(nodeNumber, position);
}

/**
* Получить номер узла в сетке по его положению
* @param position - предполагаемая позиция узла в сетке
* @param nodeNumber - номер найденного узла
* @return признак, существует ли узел с данным положением
*/
bool MeshDataProvider::GetNodeByPosition
	(
		const int position[3],
		int& nodeNumber
	)	const
{
	if (_grid == nullptr)
	{
		return false;
	}

	return _grid->GetNodeByPosition(position, nodeNumber);
}

/**
* Получить координаты узла по его положению
* @param position - предполагаемая позиция узла в сетке
* @param coordinates - координаты найденного узла
* @return признак, существует ли узел с данным положением
*/
void MeshDataProvider::GetNodeCoordinatesByPosition
	(
		const int position[3],
		double coordinates[3]
	)	const
{
	if (_grid == nullptr)
	{
		return;
	}
	_grid->GetNodeCoordinatesByPosition(position, coordinates);
}

/**
* Получить положение узла сетки по его координатам
* @param supposedCoordinates - предполагаемые координаты узла
* @param position - положение узла в сетке
*/
void MeshDataProvider::GetNodePositionByCoordinates
	(
		const double supposedCoordinates[3],
		int position[3]
	)	const
{
	if (_grid == nullptr)
	{
		return;
	}
	_grid->GetNodePositionByCoordinates
		(
			supposedCoordinates,
			position
		);
}

/**
* Определить, может ли принадлежать точка с данной координатой
* ограничивающему параллелепипеду сетки
* @param direction - индекс направления оси, по которой берется координата
* @param coordinate - координата вдоль оси
*/
bool MeshDataProvider::IsCorrectCoordinate
	(
		int direction,
		double coordinate
	)	const
{
	if (_grid == nullptr)
	{
		return false;
	}

	return _grid->IsCorrectCoordinate(direction, coordinate);
}

//нумерация элементов для решателя
int MeshDataProvider::EnumerationMesh()
{
	int index;
	int nx,ny,nz;
	index = 1;

	_grid->GetGridSize(&nx, &ny, &nz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++) 
			{
				int tmp = (*_grid)[x][y][z];

				if (tmp > 0)
				{
					(*_grid)[x][y][z] = index;
					index++;
				}
			}
		}
	}
	_nodesCount = index - 1;

	return _nodesCount;
}

/**
* Освободить память, выделенную под данные сетки
*/
void MeshDataProvider::ReleaseMeshData()
{
	delete _grid;
	for (BoundaryCondition* boundaryCondition : _surfaces)
	{
		delete boundaryCondition;
	}
	//delete _freeSolverGridParams;
	for (BoundaryNormal* i : _normals)
	{
		delete i;
	}
}

void MeshDataProvider::LoadNormals(ifstream &ifs)
{
	try
	{
		// чтение нормалей
		// кол-во нормалей
		int sizes = 0;
		ifs.read(reinterpret_cast<char*>(&sizes), sizeof(int));

		// читаем каждую нормаль
		for (int i = 0; i < sizes; i++)
		{
			BoundaryNormal *normal = new BoundaryNormal();
			unsigned long tmp_ulong;
			ifs.read(reinterpret_cast<char*>(&tmp_ulong), sizeof(tmp_ulong));
			normal->SetPointNumber(tmp_ulong);
			VertexPoint pnt;
			double tmp_double;
			ifs.read(reinterpret_cast<char*>(&tmp_double), sizeof(tmp_double));
			pnt.x = tmp_double;
			ifs.read(reinterpret_cast<char*>(&tmp_double), sizeof(tmp_double));
			pnt.y = tmp_double;
			ifs.read(reinterpret_cast<char*>(&tmp_double), sizeof(tmp_double));
			pnt.z = tmp_double;
			normal->SetUnitaryNormal(pnt);
			_normals.push_back(normal);
		}
	}
	catch (const std::exception&)
	{
		for (size_t i = 0; i < _normals.size(); i++)
		{
			delete _normals[i];
		}
		_normals.clear();
	}
}

//void MeshDataProvider::LoadFreeSolverGridParams(ifstream &ifs)
//{
	// чтение грид парпмсов
	//delete _freeSolverGridParams;
	//_freeSolverGridParams = new FreeSolverGridParams();
	//_freeSolverGridParams->LoadBinary(ifs);
//}