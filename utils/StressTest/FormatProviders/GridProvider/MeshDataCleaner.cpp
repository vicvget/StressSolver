#include "MeshDataCleaner.h"

#include "MeshDataProvider.h"
#include "MeshDataSolverFormatter.h"
#include "RLCControlWriter.h"
//#include "../../Fcore/wrappers/FileRoutines.h"

#include <algorithm>


// class MeshDataCleaner

/**
* Конструктор по умолчанию
*/
MeshDataCleaner::MeshDataCleaner()
	:
		_deletionPolicyType(DeletionPolicyTypes::DeleteAllExceptBiggest)
{

}

/**
* Деструктор
*/
MeshDataCleaner::~MeshDataCleaner()
{

}

/**
* Установить наименование файла со структурой сетки
* @param gridFileName - наименование файла со структурой сетки
*/
void MeshDataCleaner::SetGridFileName
	(
		const string& gridFileName
	)
{
	if (!gridFileName.empty())
	{
		_gridFileName = gridFileName;
	}
}

/**
* Установить тип политики удаления изолированных подобластей сетки
* @param deletionPolicyType - устанавливаемый тип политики удаления подобластей сетки
*/
void MeshDataCleaner::SetDeletionPolicyType
	(
		DeletionPolicyType deletionPolicyType
	)
{
	_deletionPolicyType = deletionPolicyType;
}

/**
* Очистить сетку от изолированных областей
* @return признак успешной (true) или неуспешной (false) очистки
*/
bool MeshDataCleaner::CleanGrid() const
{
	MeshDataProvider meshProvider(_gridFileName.c_str());

	if (!meshProvider.IsGridLoaded())
	{
		return false;
	}

	OccRectilinearGrid* grid = meshProvider.GetGrid();
	BoundaryConditions boundaryConditions = meshProvider.GetBC();
	vector<BoundaryNormal*> boundaryNormals = meshProvider.GetNormals();

	CleanGrid(_deletionPolicyType, grid, boundaryConditions, boundaryNormals);
//#define DEBUG_DIRTY
#ifdef DEBUG_DIRTY
	string dirtyGridFileName;

	dirtyGridFileName = fs::CombinePath
		(
			fs::SplitDirFromPath(_gridFileName),
			fs::SplitFileFromExt(fs::SplitFileFromPath(_gridFileName)) + "_dirty.rlc"
		);

	fs::RenameFile(_gridFileName.c_str(), dirtyGridFileName.c_str());
#endif

	RLCControlWriter meshWriter(grid, boundaryConditions);// , meshProvider.GetFreeSolverGridParams(), boundaryNormals);

	if (!meshWriter.DumpMeshToFile(_gridFileName.c_str()))
	{
		meshWriter.ReleaseMeshData();

		return false;
	}
	meshWriter.ReleaseMeshData();

	return true;
}

/**
* Очистить сетку от изолированных областей
* @param deletionPolicyType - тип политики удаления изолированных подобластей сетки
* @param grid - сетка, для которой производится очистка
* @param boundaryConditions - список граничных условий сетки, для которой производится очистка
*/
// static
void MeshDataCleaner::CleanGrid
	(
		DeletionPolicyType deletionPolicyType,
		OccRectilinearGrid* grid,
		BoundaryConditions& boundaryConditions,
		vector<BoundaryNormal*>& boundaryNormals
	)
{
	MeshDataSolverFormatter meshFormatter(grid);
	Links links;

	meshFormatter.EnumerationMesh();
	meshFormatter.ExportGrid(links);

	int subdomainsCount;
	NodesBelongingList nodesBelongingList;

	subdomainsCount = FindGridSubdomains(links, nodesBelongingList);

	SubdomainsSizesList subdomainsSizesList;

	FindSubdomainsSizes(nodesBelongingList, subdomainsCount, subdomainsSizesList);
	
	SubdomainsBelongingList subdomainsBelongingList(subdomainsCount, true);

	ApplyDeletionPolicy(deletionPolicyType, subdomainsSizesList, subdomainsBelongingList);

	subdomainsSizesList.clear();

	ConvertGrid(grid, nodesBelongingList, subdomainsBelongingList);

	nodesBelongingList.clear();
	subdomainsBelongingList.clear();

	NodesIndexes newNodesIndexes;

	newNodesIndexes.resize(meshFormatter.GetCount(), NOT_A_MEMBER);
	FindNewNodesIndexes(grid, newNodesIndexes);
	ConvertBoundaryConditions(boundaryConditions, newNodesIndexes);
	ConvertBoundaryNormals(boundaryNormals, newNodesIndexes);
}

/**
* Найти несвязанные подобласти сетки
* @param links - список ссылок на соседние узлы для всех узлов сетки
* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
* @return - количество найденных подобластей сетки
*/
// static
int MeshDataCleaner::FindGridSubdomains
	(
		const Links& links,
		NodesBelongingList& nodesBelongingList
	)
{
	int nodeIndex;
	int subdomainIndex;

	nodesBelongingList.resize(links.size(), NOT_A_MEMBER);
	for (subdomainIndex = 0; FindFreeNode(nodesBelongingList, nodeIndex); subdomainIndex++)
	{
		FillSubdomain(links, nodesBelongingList, subdomainIndex, nodeIndex);
	}

	return subdomainIndex;
}

/**
* Найти свободный (не принадлежащий ни одной подобласти) узел в сетке
* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
* @param nodeIndex - индекс найденного узла
* @return - признак, найден ли узел (true) или не найден (false)
*/
// static
bool MeshDataCleaner::FindFreeNode
	(
		NodesBelongingList& nodesBelongingList,
		int& nodeIndex
	)
{
	NodesBelongingList::const_iterator nodeIterator;

	nodeIterator = std::find(nodesBelongingList.begin(), nodesBelongingList.end(), NOT_A_MEMBER);
	if (nodeIterator == nodesBelongingList.end())
	{
		return false;
	}
	nodeIndex = (int)(nodeIterator - nodesBelongingList.begin());

	return true;
}

/**
* Заполнить очередную подобласть сетки
* @param links - список ссылок на соседние узлы для всех узлов сетки
* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
* @param subdomainIndex - номер подобласти, которую необходимо заполнить
* @param startNodeIndex - индекс узла, с которого следует начать заполнение подобласти
*/
// static
void MeshDataCleaner::FillSubdomain
	(
		const Links& links,
		NodesBelongingList& nodesBelongingList,
		int subdomainIndex,
		int startNodeIndex
	)
{
	NodesIndexes fronts[2];

	fronts[0].push_back(startNodeIndex);
	nodesBelongingList[startNodeIndex] = subdomainIndex;
	for (int frontIndex = 1; ; frontIndex++)
	{
		int currentFrontIndex = frontIndex % 2;
		int previousFrontIndex = 1 - currentFrontIndex;

		NodesIndexes& currentFront = fronts[currentFrontIndex];
		NodesIndexes& previousFront = fronts[previousFrontIndex];

		if (previousFront.empty())
		{
			break;
		}
		for (int nodeIndex : previousFront)
		{
			const NodeLinks& nodeLinks = links[nodeIndex];
			
			for (int linkIndex = 0; linkIndex < NodeLinks::LinksCount(); linkIndex++)
			{
				int neighbourNodeIndex = nodeLinks[linkIndex];

				if (neighbourNodeIndex < 0)
				{
					continue;
				}

				int& nodeBelonging = nodesBelongingList[neighbourNodeIndex];

				if (nodeBelonging >= 0)
				{
					continue;
				}
				nodeBelonging = subdomainIndex;
				currentFront.push_back(neighbourNodeIndex);
			}
		}
		previousFront.clear();
	}
}

/**
* Найти размеры подобластей сетки
* @param nodesBelongingList - список признаков принадлежности узлов той или иной подобласти сетки
* @param subdomainsCount - количество подобластей сетки
* @param subdomainsSizesList - сформированный массив размеров подобластей сетки
*/
// static
void MeshDataCleaner::FindSubdomainsSizes
	(
		const NodesBelongingList& nodesBelongingList,
		int subdomainsCount,
		SubdomainsSizesList& subdomainsSizesList
	)
{
	subdomainsSizesList.resize(subdomainsCount, 0);
	for (int subdomainIndex : nodesBelongingList)
	{
		++subdomainsSizesList[subdomainIndex];
	}
}

/**
* Применить политику удаления изолированных подобластей
* @param applyingPolicyType - применяемый тип политики удаления изолированных подобластей
* @param sizesList - список размеров изолированных подобластей сетки
* @param belongingList - список флагов принадлежности изолированных подобластей чистой сетке
*/
// static
void MeshDataCleaner::ApplyDeletionPolicy
	(
		DeletionPolicyType applyingPolicyType,
		const SubdomainsSizesList& sizesList,
		SubdomainsBelongingList& belongingList
	)
{
	switch (applyingPolicyType)
	{
	case DeletionPolicyTypes::DeleteSmallSubdomains:
		DeletionPolicy<DeletionPolicyTypes::DeleteSmallSubdomains>::Apply
			(
				sizesList,
				belongingList
			);
		break;

	case DeletionPolicyTypes::DeleteAllExceptBiggest:
		DeletionPolicy<DeletionPolicyTypes::DeleteAllExceptBiggest>::Apply
			(
				sizesList,
				belongingList
			);
		break;
	}
}

/**
* Оставить в сетке только те элементы, которые принадлежат отмеченным подобластям
* @param grid - сетка, над которой производится преобразование
* @param nodesBelongingList - массив принадлежности узлов той или иной подобласти сетки
* @param subdomainsBelongingList - список признаков принадлежности подобластей сетке
*/
// static
void MeshDataCleaner::ConvertGrid
	(
		OccRectilinearGrid* grid,
		const NodesBelongingList& nodesBelongingList,
		const SubdomainsBelongingList& subdomainsBelongingList
	)
{
	int nx;
	int ny;
	int nz;

	grid->GetGridSize(&nx, &ny, &nz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				int& nodeNumber = (*grid)[x][y][z];

				if (nodeNumber <= 0)
				{
					continue;
				}
				
				int nodeIndex = nodeNumber - 1;
				int subdomainIndex = nodesBelongingList[nodeIndex];

				if (subdomainsBelongingList[subdomainIndex])
				{
					continue;
				}
				nodeNumber = 0;
			}
		}
	}
}

/**
* Найти новые значения индексов узлов сетки
* @param grid - сетка, над которой производится преобразование
* @param newNodesIndexes - список новых индексов сетки
*/
// static
void MeshDataCleaner::FindNewNodesIndexes
	(
		const OccRectilinearGrid* grid,
		NodesIndexes& newNodesIndexes
	)
{
	int nx;
	int ny;
	int nz;
	int newNodeIndex = 0;
	
	grid->GetGridSize(&nx, &ny, &nz);
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				int nodeNumber = (*grid)[x][y][z];

				if (nodeNumber <= 0)
				{
					continue;
				}

				int nodeIndex = nodeNumber - 1;
				
				newNodesIndexes[nodeIndex] = newNodeIndex;
				newNodeIndex++;
			}
		}
	}
}

/**
* Оставить в граничном условии сетки только те элементы, которые принадлежат отмеченным подобластям
* @param boundaryCondition - граничное условие, над которым производится преобразование
* @param newNodesIndexes - список новых индексов сетки
*/
// static
void MeshDataCleaner::ConvertBoundaryCondition
	(
		BoundaryCondition*& boundaryCondition,
		const NodesIndexes& newNodesIndexes
	)
{
	BoundaryCondition* newBoundaryCondition = new BoundaryCondition();

	newBoundaryCondition->SetID(boundaryCondition->GetID());
	newBoundaryCondition->SetVar(boundaryCondition->GetVar());
	newBoundaryCondition->SetName(boundaryCondition->GetName());
	for (size_t pointIndex = 0; pointIndex < boundaryCondition->GetNumBP(); pointIndex++)
	{
		int nodeNumber = boundaryCondition->GetPoint(pointIndex);
		int nodeIndex = nodeNumber - 1;
		int newNodeIndex = nodeIndex >= 0 ? newNodesIndexes[nodeIndex] : -1;

		if (newNodeIndex < 0)
		{
			continue;
		}

		int newNodeNumber = newNodeIndex + 1;

		newBoundaryCondition->AddPoint(newNodeNumber);
	}
	delete boundaryCondition;
	boundaryCondition = newBoundaryCondition;
}

/**
* Оставить в граничных условиях сетки только те элементы, которые принадлежат отмеченным подобластям
* @param boundaryConditions - список граничных условий, над которым производится преобразование
* @param newNodesIndexes - список новых индексов сетки
*/
// static
void MeshDataCleaner::ConvertBoundaryConditions
	(
		BoundaryConditions& boundaryConditions,
		const NodesIndexes& newNodesIndexes
	)
{
	for (BoundaryCondition*& boundaryCondition : boundaryConditions)
	{
		ConvertBoundaryCondition(boundaryCondition, newNodesIndexes);
	}
}

/**
* Оставить в граничных нормалях сетки только те элементы, которые принадлежат отмеченным подобластям, а для оставленных пересчитать индексы
* @param boundaryNormals - список нормалей, над которыми производится преобразование
* @param newNodesIndexes - список новых индексов сетки
*/
void MeshDataCleaner::ConvertBoundaryNormals
	(
		vector<BoundaryNormal*>& boundaryNormals,
		const NodesIndexes& newNodesIndexes
	)
{
	vector<BoundaryNormal*> newBoundaryNormals;

	for (BoundaryNormal*& boundaryNormal : boundaryNormals)
	{
		int pointNumber = boundaryNormal->GetPointNumber();
		int pointIndex = pointNumber - 1;
		int newPointIndex = pointIndex >= 0 ? newNodesIndexes[pointIndex] : -1;

		if (newPointIndex >= 0)
		{
			BoundaryNormal* newBoundaryNormal = new BoundaryNormal();

			newBoundaryNormal->SetUnitaryNormal(boundaryNormal->GetUnitaryNormal());
			newBoundaryNormal->SetPointNumber(newPointIndex + 1);
			newBoundaryNormals.push_back(newBoundaryNormal);
		}
	}
	for (BoundaryNormal* boundaryNormal : boundaryNormals)
	{
		delete boundaryNormal;
	}
	boundaryNormals = newBoundaryNormals;
}