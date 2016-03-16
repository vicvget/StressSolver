#include "ControlReader.h"

using std::ios;


ControlReader::ControlReader(void)
{
}

ControlReader::~ControlReader(void)
{
}

int ControlReader::ReadMeshFromFile1(const char *fn)
{
	ifstream is;
	is.open (fn, ios::binary );
	if(is.is_open())
	{
		if(fillGrid(&is) == NULL)
		{
			is.close();
			return -1;
		}
		is.close();
		return 1;
	}
	return 0;
}

int ControlReader::ReadMeshFromFile1(ifstream &ifs)
{
	if (ifs.is_open())
	{
		if (fillGrid(&ifs) == NULL)
		{
			return -1;
		}
		return 1;
	}
	return 0;
}

int ControlReader::ReadMeshFromFile(const char *fn)
{
	// TODO: res инициализиурется не во всех путях ветвления
	int res;
	string rlcFileName(fn);
	if (RLCHeaderFileProvider::IsNewRLCFormat(rlcFileName))
	{
		RLCHeader* rlcHeader = RLCHeaderFileProvider::LoadRLCHeader(rlcFileName);		
		ifstream ifs(rlcFileName, std::ios_base::binary);		
		if (rlcHeader->_meshPos != 0)
		{
			ifs.seekg(rlcHeader->_meshPos);
			ReadMeshFromFile1(ifs);
		}		
		if (rlcHeader->_normalsPos != 0)
		{
			ifs.seekg(rlcHeader->_normalsPos);
			LoadNormals(ifs);
		}
		if (rlcHeader->_freeSolverGridParamsPos != 0)
		{
			ifs.seekg(rlcHeader->_freeSolverGridParamsPos);
			//LoadFreeSolverGridParams(ifs);
		}
		ifs.close();
	}
	else
	{
		res = ReadMeshFromFile1(fn);
	}
	return res;
}
