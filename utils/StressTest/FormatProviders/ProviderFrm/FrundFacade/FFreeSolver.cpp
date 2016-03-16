#include "FFreeSolver.h"

FFreeSolver::FFreeSolver()
{
	_solverType = SolverTypes::ST_Thermal;
	Init();
}

FFreeSolver::FFreeSolver(const char *fileName)
{	
	_fileId = fileName;
	Load(fileName);
}

FFreeSolver::FFreeSolver(string& name, string& fileName ,SolverTypes solverType)
	: _solverType(solverType)
{
	_name = name;
	_fileId = fileName;
	Init();
}

void FFreeSolver::Init()
{
	_isEnabled = true;
	_intParamsId = 0;
	_specialParamsId = _solverType; // первые SOLVER_TYPES_NUMBER - default
	_readIco = false;
	_writeIco = false;
}

void FFreeSolver::Load(const char* path)
{
	ifstream ifs_txt(path);
	Load(ifs_txt);
	ifs_txt.close();
}

int FFreeSolver::Load(ifstream& stream)
{
	fs::ReadLineValue(stream, _name);
	fs::ReadLineValue(stream, _number);
	int tmp_int;
	fs::ReadLineValue(stream, tmp_int);
	_solverType = (SolverTypes)tmp_int;
	fs::ReadLineValue(stream, _gridFileName);
	fs::ReadLineValue(stream, _isEnabled);
	fs::ReadLineValue(stream, _intParamsId);
	fs::ReadLineValue(stream, _specialParamsId);
	if (fs::ReadLineValueNoException(stream, _readIco))
	{
		if (fs::ReadLineValueNoException(stream, _writeIco))
			std::cout << "ERROR in MPH after read ico flag \n";
	}
	_gridParams.Load(stream);
	return 1;
}

void FFreeSolver::Save(ofstream& ofs) const
{
	fs::WriteCommentedLine(ofs, _name, "name of free solver");
	fs::WriteCommentedLine(ofs, _number, "number of free solver");
	int tmp_int = (int)_solverType;
	fs::WriteCommentedLine(ofs, tmp_int, "type of free solver");
	fs::WriteCommentedLine(ofs, _gridFileName, "grid params file name");
	int enabled = _isEnabled ? 1 : 0;
	fs::WriteCommentedLine(ofs, enabled, "flag enabled");
	fs::WriteCommentedLine(ofs, _intParamsId, "int params id");	
	fs::WriteCommentedLine(ofs, _specialParamsId, "special params id");
	fs::WriteCommentedLine(ofs, _readIco, "read from <bodyindex>_<solverId>.ics if 1");
	fs::WriteCommentedLine(ofs, _writeIco, "write to <bodyindex>_<solverId>.ics if 1");
	_gridParams.Save(ofs);
}
//
//void FFreeSolver::Output(ofstream& stream, int stage)
//{
//	int a = 0;
//}

void FFreeSolver::SaveByIndex() const
{
	ofstream ofs(_fileId);
	Save(ofs);
	ofs.close();
}