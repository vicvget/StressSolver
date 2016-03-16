#ifndef ControlReader_H
#define ControlReader_H
//#pragma GCC visibility push(hidden)
#include "MeshDataProvider.h"
#include "RLCHeader.h"
#include <iostream>
#include <fstream>

/*
Абстрактный класс для классов чтения днных из файлов
*/
class ControlReader:public MeshDataProvider
{
protected:
	virtual OccRectilinearGrid* fillGrid(ifstream *is)=0;

public:
	ControlReader(void);
	~ControlReader(void);
	/** Считывает инфу из старого формата RLC 
	* Р¤СѓРЅРєС†РёСЏ С‡РёС‚Р°РµС‚ Р±РёРЅР°СЂРЅС‹Р№ С„Р°Р№Р»
	* @param const char *fn РёРјСЏ С„Р°Р№Р»Р°
	* @return РІРѕР·РІСЂР°С‰РµС‚ 1 РµСЃР»Рё С‡С‚РµРЅРёРµ РїСЂРѕРёР·РІРµРґРµРЅРѕ СѓСЃРїРµС€РЅРѕ
	*/
	int ReadMeshFromFile1(const char *fn);
	int ReadMeshFromFile1(ifstream &ifs);

	/**
	* Р¤СѓРЅРєС†РёСЏ С‡РёС‚Р°РµС‚ Р±РёРЅР°СЂРЅС‹Р№ С„Р°Р№Р»
	* @param const char *fn РёРјСЏ С„Р°Р№Р»Р°
	* @return РІРѕР·РІСЂР°С‰РµС‚ 1 РµСЃР»Рё С‡С‚РµРЅРёРµ РїСЂРѕРёР·РІРµРґРµРЅРѕ СѓСЃРїРµС€РЅРѕ
	*/
	int ReadMeshFromFile(const char *fn);
};
//#pragma GCC visibility pop
#endif