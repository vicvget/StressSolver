#ifndef RLCControlWriter_H
#define RLCControlWriter_H
//#pragma GCC visibility push(hidden)
#include "MeshDataProvider.h"
#include "BoundaryCondition.h"
#include "RLCHeader.h"
#include <iostream>
#include <fstream>


/*Класс для записи сеток в сжатый формат*/
class RLCControlWriter:public MeshDataProvider
{
public:
	RLCControlWriter(OccRectilinearGrid *grid);
	RLCControlWriter(OccRectilinearGrid *grid,vector<BoundaryCondition *> Surfaces);
//	RLCControlWriter(OccRectilinearGrid *grid, const vector<BoundaryCondition *> &surfaces, FreeSolverGridParams *freeSolverGridParams, const vector<BoundaryNormal*> &boundaryNormals);
	~RLCControlWriter(void);
	/**
	* Функция записывает бинарный файл со сжатым форматом хранения сетки
	* @param const char *fn имя файла
	* @return возвращет 1 если запись успешна
	*/
	int DumpMeshToFile(const char *fn);

	/** Считывает инфу из старого формата RLC 
	* Функция записывает бинарный файл со сжатым форматом хранения сетки
	* @param const char *fn имя файла
	* @return возвращет 1 если запись успешна
	*/
	int DumpMeshToFile1(const char *fn);
	int DumpMeshToFile1(ofstream &ofs);
	int SaveMeshToFile(ofstream &ofs);

	/* Saves normals for boundary points
	* @fn - file with normals
	*/
	int DumpBoundaryNormalsToFile(const char *fn);
	int DumpBoundaryNormalsToFile(ofstream &ofs);

	/* Saves free solver grid params
	* @fn - file with grid params
	*/
	int DumpFreeSolverGridParamsToFile(const char *fn);
	int DumpFreeSolverGridParamsToFile(ofstream &ofs);
};
//#pragma GCC visibility pop
#endif