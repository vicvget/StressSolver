#include "FMultiphysics.h"


/** Загрузка из текстового потока
* @param ifs - входной текстовый поток
*/
void FHeatTransferParams::Load(ifstream& ifs)
{
	ifs >> time >> timeStep >> recordStep;
	fs::ReadLineString(ifs, fileGrid);
	for(int i = 0; i < _numberOfParams; i++)
	{
		ifs >> p.params[i];
	}
}

/** Сохранение в текстовый поток
* @param ofs - выходной текстовый поток
*/
void FHeatTransferParams::Save(ofstream& ofs) const
{
	ofs << time << ' ' << timeStep << ' ' << recordStep << std::endl;
	ofs << fileGrid << std::endl;
	for(int i = 0; i < _numberOfParams; i++)
	{
		ofs << p.params[i];
	}
}

/** Загрузка из текстового потока
* @param ifs - входной текстовый поток
*/
void FStressDeformationParams::Load(ifstream& ifs)
{
	ifs >> time >> timeStep >> recordStep;
	fs::ReadLineString(ifs, fileGrid);
	ifs >> K >> E;	
	// TODO: load points
}

/** Сохранение в текстовый поток
* @param ofs - выходной текстовый поток
*/
void FStressDeformationParams::Save(ofstream& ofs) const
{
	ofs << time << ' ' << timeStep << ' ' << recordStep << std::endl;
	ofs << fileGrid << std::endl;
	ofs << this->E << ' ' << this->K << std::endl;
	// TODO: save points
}

FMultiphysics::FMultiphysics
	(
		ifstream& ifs
	)
	:
		// TODO: сделать конструктор по умолчанию в классе FStressDeformationParams и убрать эту строчку
		_stressDeformationParams()
{
	Load(ifs);
}

/** Загрузка из текстового потока
* @param ifs - входной текстовый поток
*/
void FMultiphysics::Load(ifstream& ifs)
{
	int numberOfSolvers, solverType;
	ifs >> numberOfSolvers;
	_solverTypes.clear();
	for(int i = 0; i < numberOfSolvers; i++)
	{
		ifs >> solverType;
		_solverTypes.insert(solverType);
	}
	int solverCount = 0;
	while((solverCount < numberOfSolvers) && !ifs.eof())
	{
		ifs >> solverType;
		switch(solverType)
		{
		case MPHHeatTransfer:
			_heatTransferParams.Load(ifs);
			break;
		case MPHStressDeformation:
			_heatTransferParams.Load(ifs);
			break;
		}
		solverCount++;
	}

}
/** Сохранение в текстовый поток
* @param ofs - выходной текстовый поток
*/
void FMultiphysics::Save(ofstream& ofs)
{		
	set<int>::iterator it = _solverTypes.begin();
	while(it != _solverTypes.end())
	{
		switch(*it)
		{
		case MPHHeatTransfer:
			ofs << MPHHeatTransfer << std::endl;
			_heatTransferParams.Save(ofs);
			break;
		case MPHStressDeformation:
			ofs << MPHStressDeformation << std::endl;
			_stressDeformationParams.Save(ofs);
			break;
		}
	}
}
