#include "RLCControlWriter.h"


using std::ofstream;
using std::ios;


RLCControlWriter::RLCControlWriter(OccRectilinearGrid *grid):MeshDataProvider(grid)
{
	EnumerationMesh();
}


RLCControlWriter::RLCControlWriter(OccRectilinearGrid *grid,vector<BoundaryCondition *> Surfaces):MeshDataProvider(grid,Surfaces)
{
	EnumerationMesh();
}


//RLCControlWriter::RLCControlWriter(OccRectilinearGrid *grid, const vector<BoundaryCondition *> &surfaces, FreeSolverGridParams *freeSolverGridParams, const vector<BoundaryNormal*> &boundaryNormals) : MeshDataProvider(grid, surfaces)	
//{
//	_freeSolverGridParams = freeSolverGridParams;
//	_normals = boundaryNormals;
//	EnumerationMesh();
//}

RLCControlWriter::~RLCControlWriter(void)
{
}

//сохранение сетки в файл *.rlc
int RLCControlWriter::DumpMeshToFile1(const char *fn)
{
	ofstream ofs(fn , ios::out | ios::binary);
	int res = SaveMeshToFile(ofs);
	ofs.close();
	return res;
}

int RLCControlWriter::SaveMeshToFile(ofstream &ofs)
{
	if (ofs.is_open())
	{
		double xs, ys, zs;
		ofs.write(reinterpret_cast<char*>(&_nodesCount), sizeof(int));
		_grid->GetStartPoint(&xs, &ys, &zs);
		ofs.write(reinterpret_cast<char*>(&xs), sizeof(double));
		ofs.write(reinterpret_cast<char*>(&ys), sizeof(double));
		ofs.write(reinterpret_cast<char*>(&zs), sizeof(double));
		//запись точек
		int nx, ny, nz;
		_grid->GetGridSize(&nx, &ny, &nz);

		//запись заголовка
		//запись количества ячеек по осям
		ofs.write(reinterpret_cast<char*>(&nx), sizeof(int));
		ofs.write(reinterpret_cast<char*>(&ny), sizeof(int));
		ofs.write(reinterpret_cast<char*>(&nz), sizeof(int));
		//запись размера ячеек
		ofs.write(reinterpret_cast<char*>(&_dx), sizeof(double));
		ofs.write(reinterpret_cast<char*>(&_dy), sizeof(double));
		ofs.write(reinterpret_cast<char*>(&_dz), sizeof(double));

		//порядок обхода z y x
		for (int z = 0; z<nz; z++)
		{
			for (int y = 0; y<ny; y++)
			{
				for (int x = 0; x<nx; x++)
				{
					bool tmp = ((*_grid)[x][y][z]>0);
					ofs.write(reinterpret_cast<char*>(&tmp), sizeof(bool));
				}
			}
		}
		//запись граничных условий
		//const char *buffer;
		ofstream textofs;
#ifdef DEBUG_RLC
		textofs.open("debug");
#endif
		int sizes = static_cast<int>(_surfaces.size());
		ofs.write(reinterpret_cast<char*>(&sizes), sizeof(int));
		if (textofs.is_open())
			textofs << sizes << std::endl;
		for (int i = 0; i<sizes; i++)
		{
			BoundaryCondition *face = _surfaces.at(i);

			string id = face->GetID();
			int tmp = static_cast<int>(id.length());
			ofs.write((char *)&tmp, sizeof(int));
			ofs.write(id.c_str(), tmp);
			if (textofs.is_open())
				textofs << id << std::endl;

			string var = face->GetVar();
			tmp = static_cast<int>(var.length());
			ofs.write((char *)&tmp, sizeof(int));
			ofs.write(var.c_str(), tmp);
			if (textofs.is_open())
				textofs << var << std::endl;

			string name = face->GetName();
			tmp = static_cast<int>(name.length());
			ofs.write((char *)&tmp, sizeof(int));
			ofs.write(name.c_str(), tmp);
			if (textofs.is_open())
				textofs << name << std::endl;

			int numBP = static_cast<int>(face->GetNumBP());
			ofs.write(reinterpret_cast<char*>(&numBP), sizeof(int));
			if (textofs.is_open())
				textofs << numBP << std::endl;

			for (int j = 0; j<numBP; j++)
			{
				int pnt = face->GetPoint(j);
				ofs.write(reinterpret_cast<char*>(&pnt), sizeof(int));
				if (textofs.is_open())
					textofs << pnt << std::endl;
			}
		}
		if (textofs.is_open())
			textofs.close();
		//_ofs.close();
		return 1;
	}
	return 0;
}

int RLCControlWriter::DumpMeshToFile(const char *fn)
{
	RLCHeader header;
	header._version = -1;
	ofstream ofs(fn, ios::binary);
	RLCHeaderFileProvider::SaveRLCHeader(ofs, &header);
	if ((int)ofs.tellp() >= 0)
	{
		header._meshPos = static_cast<unsigned int>(ofs.tellp());
	}
	SaveMeshToFile(ofs);
	if (_normals.size() != 0)
	{
		if ((int)ofs.tellp() >= 0)
		{
			header._normalsPos = static_cast<unsigned int>(ofs.tellp());
		}
		DumpBoundaryNormalsToFile(ofs);
	}
	//if (_freeSolverGridParams != nullptr)
	//{
	//	if ((int)_ofs.tellp() >= 0)
	//	{
	//		header._freeSolverGridParamsPos = static_cast<unsigned int>(_ofs.tellp());
	//	}
	//	DumpFreeSolverGridParamsToFile(_ofs);
	//}
	ofs.seekp(ios::beg);
	RLCHeaderFileProvider::SaveRLCHeader(ofs, &header);
	ofs.close();
	return 1;
}

int RLCControlWriter::DumpBoundaryNormalsToFile(const char *fn)
{
	ofstream ofs(fn, ios::binary);
	DumpBoundaryNormalsToFile(ofs);
	ofs.close();
	return 1;
}

int RLCControlWriter::DumpBoundaryNormalsToFile(ofstream &ofs)
{
	// считаем и записываем кол-во граничных точек
	int pointsCount = static_cast<int>(_normals.size());
	ofs.write(reinterpret_cast<char*>(&pointsCount), sizeof(pointsCount));
	// пишем нормали для каждой точки
	for (int i = 0; i < pointsCount; i++)
	{
		// норомаль представляет собой номер точки, к которой она относится, и координаты единичного вектора вдоль OX, OY и OZ
		unsigned long tmp_ulong = _normals.at(i)->GetPointNumber();
		ofs.write(reinterpret_cast<char*>(&tmp_ulong), sizeof(tmp_ulong));
		double pnt = _normals.at(i)->GetUnitaryNormal().x;
		ofs.write(reinterpret_cast<char*>(&pnt), sizeof(pnt));
		pnt = _normals.at(i)->GetUnitaryNormal().y;
		ofs.write(reinterpret_cast<char*>(&pnt), sizeof(pnt));
		pnt = _normals.at(i)->GetUnitaryNormal().z;
		ofs.write(reinterpret_cast<char*>(&pnt), sizeof(pnt));
	}
	return 1;
}

int RLCControlWriter::DumpFreeSolverGridParamsToFile(const char *fn)
{
	ofstream ofs(fn, ios::binary);
	//_freeSolverGridParams->SaveBinary(_ofs);
	ofs.close();
	return 1;
}

int RLCControlWriter::DumpFreeSolverGridParamsToFile(ofstream &ofs)
{	
	//_freeSolverGridParams->SaveBinary(_ofs);	
	return 1;
}