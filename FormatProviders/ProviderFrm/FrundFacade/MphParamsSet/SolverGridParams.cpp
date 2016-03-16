#include "SolverGridParams.h"

#include "../../../GridProvider/MeshDataSolverFormatter.h"
#include "../../../../fcore/wrappers/FileRoutines.h"


using std::fixed;


SolverGridParams::SolverGridParams()
{
	Init();
}

SolverGridParams::SolverGridParams(const string& name)
	:
		ParamsSetBase(name)
{
	Init();
}

void SolverGridParams::Init()
{
	//_pathCAD = "";
	//_pathGrid = "";
	_gridStep = 1;
	_bcMapperId = 0;
	_isMirrored = false;
}

void SolverGridParams::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	int isMirrored = _isMirrored ? 1 : 0;
	// точность вывода шага сетки (в знаках после запятой)
	// TODO: выставить высокую точность для других double-параметров специальных решателей
	ofs.precision(10);
	//ofs	<< _pathCAD << std::endl
	//	<< _pathGrid << std::endl
	//	<< fixed << _gridStep << std::endl // для вывода шага сетки нужна большая точность!!!
	//	<< _bcMapperId << std::endl
	//	<< isMirrored << std::endl
	//	<< CountTransformPoints() << std::endl;

	ofs	<< _pathCAD << std::endl
		<< _pathGrid << std::endl
		<< fixed << _gridStep << std::endl;

	fs::WriteCommentedLine(ofs, _bcMapperId, "bc mapper id");
	fs::WriteCommentedLine(ofs, isMirrored, "flag is mirrored");
	fs::WriteCommentedLine(ofs, CountTransformPoints(), "count points");
	//	<< fixed << _gridStep << std::endl // для вывода шага сетки нужна большая точность!!!
	//	<< _bcMapperId << std::endl
	//	<< isMirrored << std::endl
	//	<< CountTransformPoints() << std::endl;

	PointsMapper::const_iterator it = _transformPointsMapper.begin();

	while (it != _transformPointsMapper.end())
	{
		ofs << it->first << std::endl;
		for (int i = 0; i < 3; i++)
		{
			ofs << it->second.coords[i] << std::endl;
		}
		++it;
	}
		
	fs::WriteCommentedLine(ofs, CountBcMappers(), "count bc mappers");
	//ofs	<< CountBcMappers() << std::endl;
	for(int i = 0; i < CountBcMappers(); i++)
		_bcMappers[i].Save(ofs);
}

void SolverGridParams::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	fs::ReadLineString(ifs,_pathCAD);
	fs::ReadLineString(ifs,_pathGrid);
	fs::ReadLineValue(ifs,_gridStep);
	fs::ReadLineValue(ifs,_bcMapperId);
	fs::ReadLineValue(ifs,_isMirrored);
	size_t pointsCount;
	
	fs::ReadLineValue(ifs,pointsCount);
	_transformPointsMapper.clear();
	
	for(int i = 0; i < pointsCount; i++)
	{
		int id;
		FGeometryPoint point;
		fs::ReadLineValue(ifs,id);
		fs::ReadLineValue(ifs,point.x);
		fs::ReadLineValue(ifs,point.y);
		fs::ReadLineValue(ifs,point.z);

		AddPoint(id, point);
	}

	size_t bcMappersCount;
	fs::ReadLineValue(ifs,bcMappersCount);
	_bcMappers.resize(bcMappersCount);
	for(int i = 0; i < CountBcMappers(); i++)	
		_bcMappers[i].Load(ifs);
}

BcMapper& SolverGridParams::CreateNewBcMapper()
{
	BcMapper mapper;
	_bcMapperId = _bcMappers.size();
	_bcMappers.push_back(mapper);
	return _bcMappers[_bcMappers.size()-1];	
}

void SolverGridParams::GetNodesCountsInBcs(vector<size_t> &fullBoundaryNodesCounts)
{
	MeshDataSolverFormatter *gridProvider = new MeshDataSolverFormatter(PathGrid().c_str());
	if (!gridProvider->IsGridLoaded())
		exceptions::ThrowMessage("Grid is not loaded");

	BcMapper &bcMapper = GetCurrentMapper();
	fullBoundaryNodesCounts.resize(bcMapper.CountBcs());
	// добавляем граничные условия для заданых поверхностей
	// добавляем граничные условия для заданых поверхностей
	for (size_t bcId = 0; bcId < bcMapper.CountBcs(); bcId++)
	{
		BcSurface &bcSurface = bcMapper.GetBc(bcId);
		if (!bcSurface.IsEnabled())
			continue;
		int *boundaryNodesIndices;
		size_t numberOfBoundaryNodes;
		gridProvider->ExportBoundary(boundaryNodesIndices, numberOfBoundaryNodes, bcSurface.GetSurfaceIds());
		fullBoundaryNodesCounts[bcId] = numberOfBoundaryNodes;
	}
}