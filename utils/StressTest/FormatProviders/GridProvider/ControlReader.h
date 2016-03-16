#ifndef ControlReader_H
#define ControlReader_H
//#pragma GCC visibility push(hidden)
#include "MeshDataProvider.h"
#include "RLCHeader.h"
#include <iostream>
#include <fstream>

/*
����������� ����� ��� ������� ������ ����� �� ������
*/
class ControlReader:public MeshDataProvider
{
protected:
	virtual OccRectilinearGrid* fillGrid(ifstream *is)=0;

public:
	ControlReader(void);
	~ControlReader(void);
	/** ��������� ���� �� ������� ������� RLC 
	* Функция читает бинарный файл
	* @param const char *fn имя файла
	* @return возвращет 1 если чтение произведено успешно
	*/
	int ReadMeshFromFile1(const char *fn);
	int ReadMeshFromFile1(ifstream &ifs);

	/**
	* Функция читает бинарный файл
	* @param const char *fn имя файла
	* @return возвращет 1 если чтение произведено успешно
	*/
	int ReadMeshFromFile(const char *fn);
};
//#pragma GCC visibility pop
#endif