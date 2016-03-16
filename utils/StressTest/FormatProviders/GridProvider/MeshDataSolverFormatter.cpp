#include "MeshDataSolverFormatter.h"

#include <algorithm>

using std::find;


// class MeshDataSolverFormatter

/**
* Конструктор
* @param gridFileName - наименование файла со структурой сетки
*/
MeshDataSolverFormatter::MeshDataSolverFormatter
	(
		const string& gridFileName
	)
	:
		MeshDataProvider(gridFileName)
{

}

/**
* Конструктор
* @param grid - структура регулярной ортогональной сетки
*/
MeshDataSolverFormatter::MeshDataSolverFormatter
	(
		OccRectilinearGrid* grid
	)
	:
		MeshDataProvider(grid)
{

}

/**
* Получить шаг сетки
* @return шаг сетки
*/
double MeshDataSolverFormatter::GetStep() const
{
	return _dx;
}

/**
* Формирует список ссылок на соседние узлы для всех узлов сетки
* @param links - сформированный список ссылок на соседние узлы для всех узлов сетки
*/
void MeshDataSolverFormatter::ExportGrid
	(
		Links& links
	)	const
{
	int sizes[3];
	int& nx = sizes[0];
	int& ny = sizes[1];
	int& nz = sizes[2];
	int coordinates[3];
	int& x = coordinates[0];
	int& y = coordinates[1];
	int& z = coordinates[2];
	int dCoordinates[3] = {0};
	const int& dx = dCoordinates[0];
	const int& dy = dCoordinates[1];
	const int& dz = dCoordinates[2];

	_grid->GetGridSize(&nx, &ny, &nz);
	links.resize(_nodesCount);
	for (z = 0; z < nz; z++)
	{
		for (y = 0; y < ny; y++)
		{
			for (x = 0; x < nx; x++)
			{
				int nodeNumber = (*_grid)[x][y][z];
				
				if (nodeNumber <= 0)
				{
					continue;
				}

				NodeLinks& nodeLinks = links[nodeNumber - 1];

				for (int direction = 0; direction < NodeLinks::DirectionsCount; direction++)
				{
					LinkDirection& linkDirection = nodeLinks[(Direction)direction];
					const int& size = sizes[direction];
					const int& coordinate = coordinates[direction];
					int& dCoordinate = dCoordinates[direction];

					for (int variant = 0; variant < LinkDirection::VariantsCount; variant++)
					{
						switch (variant)
						{
						case LinkDirection::MinusVariant:
							if (coordinate <= 0)
							{
								continue;
							}
							dCoordinate = -1;
							break;

						case LinkDirection::PlusVariant:
							if (coordinate >= size - 1)
							{
								continue;
							}
							dCoordinate = 1;
							break;
						}
						/*if
							(
								((variant == LinkDirection::MinusVariant) && (coordinate > 0)) ||
								((variant == LinkDirection::PlusVariant) && (coordinate < size - 1))
							)*/
						//dCoordinate = variant ? 1 : -1;

						int neighbourNumber = (*_grid)[x + dx][y + dy][z + dz];

						if (neighbourNumber > 0)
						{
							linkDirection[(DirectionVariant)variant] = neighbourNumber - 1;
						}
						dCoordinate = 0;
					}
				}
			}
		}
	}
}

void MeshDataSolverFormatter::ExportGrid
	(
		int*& mlink,
		double*& uzk
	)
{
	int nx;
	int ny;
	int nz;
	int kuzId = 0;
	int mlinkId = 0;
	double xs;
	double ys;
	double zs;
	EnumerationMesh();
	_grid->GetGridSize(&nx, &ny, &nz);
	_grid->GetStartPoint(&xs, &ys, &zs);
	mlink = new int[_nodesCount * 6];
	uzk = new double[_nodesCount * 3];
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				int tmp = (*_grid)[x][y][z];
				if (tmp > 0)
				{
					uzk[kuzId++] = (xs + x * _dx) / 1000;
					uzk[kuzId++] = (ys + y * _dy) / 1000;
					uzk[kuzId++] = (zs + z * _dz) / 1000;
					//if (z > 0)
					//	mlink[mlinkId++] = ((*grid)[x][y][z - 1] > 0) ? (*grid)[x][y][z - 1] : 0;
					//else
					//	mlink[mlinkId++]=0;
					//if (z < nz-1)
					//	mlink[mlinkId++] = ((*grid)[x][y][z + 1] > 0) ? (*grid)[x][y][z + 1] : 0;
					//else
					//	mlink[mlinkId++]=0;
					//if (y > 0)
					//	mlink[mlinkId++] = ((*grid)[x][y - 1][z] > 0) ? (*grid)[x][y - 1][z] : 0;
					//else
					//	mlink[mlinkId++]=0;
					//if (y < ny-1)
					//	mlink[mlinkId++] = ((*grid)[x][y + 1][z] > 0) ? (*grid)[x][y + 1][z] : 0;
					//else
					//	mlink[mlinkId++]=0;
					//if (x > 0)
					//	mlink[mlinkId++] = ((*grid)[x - 1][y][z] > 0) ? (*grid)[x - 1][y][z] : 0;
					//else
					//	mlink[mlinkId++]=0;
					//if (x < nx-1)
					//	mlink[mlinkId++] = ((*grid)[x + 1][y][z] > 0) ? (*grid)[x + 1][y][z] : 0;
					//else
					//	mlink[mlinkId++]=0;

					if (x > 0)
						mlink[mlinkId++] = ((*_grid)[x - 1][y][z] > 0) ? (*_grid)[x - 1][y][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (y > 0)
						mlink[mlinkId++] = ((*_grid)[x][y - 1][z] > 0) ? (*_grid)[x][y - 1][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (z > 0)
						mlink[mlinkId++] = ((*_grid)[x][y][z - 1] > 0) ? (*_grid)[x][y][z - 1] : 0;
					else
						mlink[mlinkId++]=0;
					if (x < nx-1)
						mlink[mlinkId++] = ((*_grid)[x + 1][y][z] > 0) ? (*_grid)[x + 1][y][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (y < ny-1)
						mlink[mlinkId++] = ((*_grid)[x][y + 1][z] > 0) ? (*_grid)[x][y + 1][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (z < nz-1)
						mlink[mlinkId++] = ((*_grid)[x][y][z + 1] > 0) ? (*_grid)[x][y][z + 1] : 0;
					else
						mlink[mlinkId++]=0;
				}			
			}
		}
	}

}

void MeshDataSolverFormatter::ExportGrid
	(
		int*& mlink
	)
{
	int nx;
	int ny;
	int nz;
	int mlinkId = 0;
	double xs;
	double ys;
	double zs;
	EnumerationMesh();
	_grid->GetGridSize(&nx, &ny, &nz);
	_grid->GetStartPoint(&xs, &ys, &zs);
	mlink = new int[_nodesCount * 6];
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				int tmp = (*_grid)[x][y][z];
				if (tmp > 0)
				{
					if (z > 0)
						mlink[mlinkId++] = ((*_grid)[x][y][z - 1] > 0) ? (*_grid)[x][y][z - 1] : 0;
					else
						mlink[mlinkId++]=0;
					if (z < nz-1)
						mlink[mlinkId++] = ((*_grid)[x][y][z + 1] > 0) ? (*_grid)[x][y][z + 1] : 0;
					else
						mlink[mlinkId++]=0;
					if (y > 0)
						mlink[mlinkId++] = ((*_grid)[x][y - 1][z] > 0) ? (*_grid)[x][y - 1][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (y < ny-1)
						mlink[mlinkId++] = ((*_grid)[x][y + 1][z] > 0) ? (*_grid)[x][y + 1][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (x > 0)
						mlink[mlinkId++] = ((*_grid)[x - 1][y][z] > 0) ? (*_grid)[x - 1][y][z] : 0;
					else
						mlink[mlinkId++]=0;
					if (x < nx-1)
						mlink[mlinkId++] = ((*_grid)[x + 1][y][z] > 0) ? (*_grid)[x + 1][y][z] : 0;
					else
						mlink[mlinkId++]=0;
				}			
			}
		}
	}
}

void MeshDataSolverFormatter::ExportBoundary
	(
		int*& mbon,
		size_t& kbon
	)	const
{
	kbon = 0;
	for (BoundaryCondition* face : _surfaces)
	{
		kbon += face->GetNumBP();
	}
	mbon = new int[kbon];
	
	kbon = 0;
	for (BoundaryCondition* face : _surfaces)
	{
		if (face->GetNumBP() > 0)
		{
			for (size_t j = 0; j < face->GetNumBP(); j++)
			{
				mbon[kbon++] = face->GetPoint(j);
			}
		}
	}
}

void MeshDataSolverFormatter::ExportBoundary
	(
		int*& mbon,
		size_t& kbon,
		const vector<string>& surfaceIds
	)	const
{
	kbon = 0;
	for (BoundaryCondition* face : _surfaces)
	{
		string name = face->GetName();
		if(std::find(surfaceIds.begin(),surfaceIds.end(),name) != surfaceIds.end())
		{
			kbon += face->GetNumBP();
		}
	}
	mbon = new int[kbon];
	
	kbon = 0;
	for (BoundaryCondition* face : _surfaces)
	{
		string name = face->GetName();
		if (std::find(surfaceIds.begin(), surfaceIds.end(), name) != surfaceIds.end())
		{
			if (face->GetNumBP() > 0)
			{
				for (size_t j = 0; j < face->GetNumBP(); j++)
				{
					mbon[kbon++] = face->GetPoint(j);
				}
			}
		}
	}
}

void MeshDataSolverFormatter::ExportBoundary
	(
		const vector<string>& surfaceIds,
		vector<int>& boundaryNodesIndices
	)	const
{
	int* boundaryNodesIndicesTmp; // индексы граничных точек
	size_t numberOfBoundaryNodes; // количество граничных точек

	// извлекаем точки из сетки
	ExportBoundary
		(
			boundaryNodesIndicesTmp,
			numberOfBoundaryNodes,
			surfaceIds
		);

	if (numberOfBoundaryNodes > 0)
	{
		boundaryNodesIndices.clear();
		boundaryNodesIndices.resize(numberOfBoundaryNodes);
		for (size_t nodeIndex = 0; nodeIndex < boundaryNodesIndices.size(); nodeIndex++)
		{
			boundaryNodesIndices.at(nodeIndex) = boundaryNodesIndicesTmp[nodeIndex] - 1;
		}
		//memcpy(&boundaryNodesIndices[0], boundaryNodesIndicesTmp, numberOfBoundaryNodes * sizeof(int));
	}
}