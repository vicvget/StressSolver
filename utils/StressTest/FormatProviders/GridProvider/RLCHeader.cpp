#include "RLCHeader.h"

//#include "../../Fcore/Exceptions/fcExceptions.h"


RLCHeader* RLCHeaderFileProvider::LoadRLCHeader(ifstream &ifs)
{
	RLCHeader* rlcHeader = new RLCHeader();		
	ifs.read(reinterpret_cast<char*>(&(rlcHeader->_version)), sizeof(rlcHeader->_version));	
	ifs.read(reinterpret_cast<char*>(&(rlcHeader->_meshPos)), sizeof(rlcHeader->_meshPos));	
	ifs.read(reinterpret_cast<char*>(&(rlcHeader->_normalsPos)), sizeof(rlcHeader->_normalsPos));	
	ifs.read(reinterpret_cast<char*>(&(rlcHeader->_freeSolverGridParamsPos)), sizeof(rlcHeader->_freeSolverGridParamsPos));	
	return rlcHeader;
}

RLCHeader* RLCHeaderFileProvider::LoadRLCHeader(const string &rlcFileName)
{
	ifstream ifs(rlcFileName, ifstream::binary);
	RLCHeader *rlcHeader = LoadRLCHeader(ifs);
	ifs.close();
	return rlcHeader;
}

void RLCHeaderFileProvider::SaveRLCHeader(ofstream &ofs, RLCHeader* rlcHeader)
{
	ofs.write(reinterpret_cast<char*>(&(rlcHeader->_version)), sizeof(rlcHeader->_version));
	ofs.write(reinterpret_cast<char*>(&(rlcHeader->_meshPos)), sizeof(rlcHeader->_meshPos));
	ofs.write(reinterpret_cast<char*>(&(rlcHeader->_normalsPos)), sizeof(rlcHeader->_normalsPos));
	ofs.write(reinterpret_cast<char*>(&(rlcHeader->_freeSolverGridParamsPos)), sizeof(rlcHeader->_freeSolverGridParamsPos));
}

bool RLCHeaderFileProvider::IsNewRLCFormat(const string &rlcFileName)
{
	ifstream ifs(rlcFileName, ifstream::binary);

	if (!ifs.is_open())
	{
		throw "FileNotFound";
		//exceptions::ThrowFileNotFound(rlcFileName);
	}

	RLCHeader* rlcHeader = LoadRLCHeader(ifs);

	return rlcHeader->_version < 0;
}
