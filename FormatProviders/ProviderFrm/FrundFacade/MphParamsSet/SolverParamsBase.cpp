#include "SolverParamsBase.h"

#include "../../../../fcore/wrappers/FileRoutines.h"

#include <iostream>


using std::cout;
using std::endl;


SolverParamsBase::SolverParamsBase(string& name, SolverTypes solverType)
: ParamsSetBase(name), _solverType(solverType)
{
	Init();
}

SolverParamsBase::SolverParamsBase(SolverTypes solverType)
: ParamsSetBase(), _solverType(solverType)
{
	Init();
}

void SolverParamsBase::Init()
{
	_isEnabled = true;
	_intParamsId = 0;
	_gridParamsId = 0;
	_specialParamsId = _solverType; // первые SOLVER_TYPES_NUMBER - default
	_solidTransparency = 0;
	_readIco = false;
	_writeIco = false;
}

void SolverParamsBase::Save(ofstream& ofs) const
{
	DefaultSave(ofs);	
	int enabled = _isEnabled ? 1 : 0;
	fs::WriteCommentedLine(ofs, enabled, "flag enabled");
	fs::WriteCommentedLine(ofs, _intParamsId, "int params id");
	fs::WriteCommentedLine(ofs, _gridParamsId, "grid params id");
	fs::WriteCommentedLine(ofs, _specialParamsId, "special params id");
	fs::WriteCommentedLine(ofs, _readIco, "read from <bodyindex>_<solverId>.ics if 1");
	fs::WriteCommentedLine(ofs, _writeIco, "write to <bodyindex>_<solverId>.ics if 1");
}

void SolverParamsBase::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	fs::ReadLineValue(ifs,_isEnabled);
	fs::ReadLineValue(ifs,_intParamsId);
	fs::ReadLineValue(ifs,_gridParamsId);
	fs::ReadLineValue(ifs,_specialParamsId);
	if(fs::ReadLineValueNoException(ifs,_readIco))
	{
		if(fs::ReadLineValueNoException(ifs,_writeIco))
			cout << "ERROR in MPH after read ico flag \n";
	}
	else
	{
		cout << "WARNING in MPH: Old format, resave will fix it only if there are no bodies with multiple solvers \n";
	}
}
