#include "FMph.h"

#include "../../../fcore/Exceptions/fcExceptions.h"

#include <vector>


using std::vector;


/** Добавляет и возвращает новый решатель
* @param type - тип решателя
* @return решатель
*/
SolverParamsBase* FMph::AddSolver(SolverTypes type)
{	
	// Тут можно создавать унаследованные решатели в зависимости от типа, 
	// поэтому тип и вынесен за пределы сохранения самого решателя
	SolverParamsBase* solver = new SolverParamsBase(type);
	_solverParams.push_back(solver);
	return _solverParams[_solverParams.size()-1];
}

/** Удаляет решатель по номеру
* @param solverId - неомер решателя
*/
void FMph::RemoveSolver(size_t solverId)
{
	if(solverId < _solverParams.size())
		_solverParams.erase(_solverParams.begin()+solverId);
	else
		exceptions::ThrowCollectionOutOfBounds("_solverParams");
}

/** Возвращает решателья по номеру
* @param solverId - неомер решателя
* @return решатель
*/
SolverParamsBase* FMph::GetSolver(size_t solverId)
{
	if(solverId > _solverParams.size())
		exceptions::ThrowCollectionOutOfBounds("_solverParams");
	return _solverParams[solverId];
}


void FMph::Save(ofstream& ofs) const
{
	fs::WriteCommentedLine(ofs, CountSolvers(), "count solvers");
	//ofs << CountSolvers() << std::endl;
	for(size_t i = 0; i < CountSolvers(); i++)
	{
		fs::WriteCommentedLine(ofs, _solverParams[i]->GetSolverType(), "solver type");
		//ofs << _solverParams[i]->GetSolverType() << std::endl;
		_solverParams[i]->Save(ofs);
	}
}

bool FMph::Load(ifstream& ifs) 
{
	int numberOfSolvers;

	fs::ReadLineValue(ifs,numberOfSolvers);

	int solverType;

	for (size_t i = 0; i < numberOfSolvers; i++)
	{
		fs::ReadLineValue(ifs,solverType);
		AddSolver((SolverTypes)solverType)->Load(ifs);
	}

	return false;
}

FMph::~FMph()
{
	for(int i = 0; i < _solverParams.size(); i++)
	{
		delete _solverParams[i];
	}
}
