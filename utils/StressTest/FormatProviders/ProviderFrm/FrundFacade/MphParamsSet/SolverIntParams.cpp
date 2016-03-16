#include "SolverIntParams.h"
#include "../../../../fcore/wrappers/FileRoutines.h"

SolverIntParams::SolverIntParams()
: ParamsSetBase()
{
	Init();
}

SolverIntParams::SolverIntParams(string &name)
: ParamsSetBase(name)
{
	Init();
}

void SolverIntParams::Init()
{
	_timeStep = 1e-3; 
	_numberOfFrames = 10;
}

void SolverIntParams::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	ofs	<< _timeStep << std::endl
		<< _numberOfFrames << std::endl;
}

void SolverIntParams::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	fs::ReadLineValue(ifs,_timeStep);
	fs::ReadLineValue(ifs,_numberOfFrames);
}
