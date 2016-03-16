#include "FMphHeader.h"

#include "../../../fcore/wrappers/FileRoutines.h"

#include <fstream>
#include <string>


using std::endl;


FMphHeader::FMphHeader():FElement()
{
	Init();
}

FMphHeader::FMphHeader(const char* path)
{
	Init();
	FElement::Load(path);
}

void FMphHeader::Init()
{
	_intParamsSets.clear();
	_gridParamsSets.clear();
	_specialParamsSets.clear();
	_bcParamsSets.clear();

	AddNewIntParams();
	AddNewGridParams();
	// TODO: до SOLVER_TYPES_NUMBER
	AddNewSpecialParams(ST_Thermal); // default = thermal
	AddNewSpecialParams(ST_StressStrain); // default = stressstrain
	// TODO: до BC_TYPES_NUMBER
	AddNewBcParams(BCT_ThermalBoundary1);
	AddNewBcParams(BCT_ThermalBoundary2);
	AddNewBcParams(BCT_ThermalBoundary3);
	AddNewBcParams(BCT_StressStrainBoundaryForce);
	AddNewBcParams(BCT_StressStrainBoundarySealing);
}

SolverSpecialParams& FMphHeader::AddNewSpecialParams(SolverTypes type)
{
	SolverSpecialParams specialParams(type);
	_specialParamsSets.push_back(specialParams);	
	return _specialParamsSets[_specialParamsSets.size()-1];
}

SolverIntParams& FMphHeader::AddNewIntParams()
{
	SolverIntParams intParams;
	_intParamsSets.push_back(intParams);	
	return _intParamsSets[_intParamsSets.size()-1];
}

SolverIntParams& FMphHeader::GetIntParams(size_t id)
{
	if(!(id < _intParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_intParamsSets");
	return _intParamsSets[id];
}

SolverGridParams& FMphHeader::GetGridParams(size_t id)
{
	if(!(id < _gridParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_gridParamsSets");
	return _gridParamsSets[id];
}

SolverSpecialParams& FMphHeader::GetSpecialParams(size_t id)
{
	if(!(id < _specialParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_specialParamsSets");
	return _specialParamsSets[id];
}

BcParams& FMphHeader::GetBcParams(size_t id)
{
	if(!(id < _bcParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_bcParams");
	return _bcParamsSets[id];
}
// const
const SolverIntParams& FMphHeader::GetIntParams(size_t id) const
{
	if(!(id < _intParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_intParamsSets");
	return _intParamsSets[id];
}

const SolverGridParams& FMphHeader::GetGridParams(size_t id) const
{
	if(!(id < _gridParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_gridParamsSets");
	return _gridParamsSets[id];
}

const SolverSpecialParams& FMphHeader::GetSpecialParams(size_t id) const
{
	if(!(id < _specialParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_specialParamsSets");
	return _specialParamsSets[id];
}

const BcParams& FMphHeader::GetBcParams(size_t id) const
{
	if(!(id < _bcParamsSets.size()))
		exceptions::ThrowCollectionOutOfBounds("_bcParams");
	return _bcParamsSets[id];
}

SolverGridParams& FMphHeader::AddNewGridParams()
{
	SolverGridParams gridParams;
	_gridParamsSets.push_back(gridParams);
	return _gridParamsSets[_gridParamsSets.size()-1];
}

BcParams& FMphHeader::AddNewBcParams(BcTypes type)
{
	BcParams bcParams(type);
	_bcParamsSets.push_back(bcParams);
	return _bcParamsSets[_bcParamsSets.size()-1];
}

// virtual
bool FMphHeader::Load(ifstream& ifs) // override
{
	Init();

	size_t tmpCount;
	size_t i;

	fs::ReadLineValue(ifs, tmpCount);
	tmpCount++;
	_intParamsSets.resize(tmpCount);
	for (i = 1; i < tmpCount; i++)
	{
		_intParamsSets[i].Load(ifs);
	}

	fs::ReadLineValue(ifs, tmpCount);
	tmpCount++;
	_gridParamsSets.resize(tmpCount);
	for (i = 1; i < tmpCount; i++)
	{
		_gridParamsSets[i].Load(ifs);
	}
	fs::ReadLineValue(ifs, tmpCount);
	tmpCount+=SOLVER_TYPES_NUMBER;
	
	SolverSpecialParams tmpParam(ST_Thermal);

	for (i = SOLVER_TYPES_NUMBER; i < tmpCount; i++)
	{
		tmpParam.Load(ifs);
		_specialParamsSets.push_back(tmpParam);
	}
	fs::ReadLineValue(ifs, tmpCount);
	tmpCount += BC_TYPES_NUMBER; // одно значение default

	// не важно какого типа создастся объект, тип потом прочитается и перезапишется
	BcParams tmpBcParams(BCT_ThermalBoundary1);

	for (i = BC_TYPES_NUMBER; i < tmpCount; i++)
	{
		tmpBcParams.Load(ifs);
		_bcParamsSets.push_back(tmpBcParams);
	}

	return true;
}

// virtual
void FMphHeader::Save(ofstream& ofs) const // override
{
	size_t i;
	//ofs << CountIntParamsSets()-1 << endl;
	fs::WriteCommentedLine(ofs, CountIntParamsSets()-1, "int params count");
	for(i = 1; i < CountIntParamsSets(); i++)
	{
		_intParamsSets[i].Save(ofs);
	}

	//ofs << CountGridParamsSets()-1 << endl;
	fs::WriteCommentedLine(ofs, CountGridParamsSets()-1, "grid params count");
	for(i = 1; i < CountGridParamsSets(); i++)
	{
		_gridParamsSets[i].Save(ofs);
	}
	
	//ofs << CountSpecialParamsSets()-SOLVER_TYPES_NUMBER << endl;
	fs::WriteCommentedLine(ofs, CountSpecialParamsSets()-SOLVER_TYPES_NUMBER, "special params count");
	for(i = SOLVER_TYPES_NUMBER; i < CountSpecialParamsSets(); i++)
	{
		_specialParamsSets[i].Save(ofs);
	}

	//ofs << CountBcParamsSets()-BC_TYPES_NUMBER << endl;
	fs::WriteCommentedLine(ofs, CountBcParamsSets()-BC_TYPES_NUMBER, "bc params count");
	for(i = BC_TYPES_NUMBER; i < CountBcParamsSets(); i++)
	{
		_bcParamsSets[i].Save(ofs);
	}
}