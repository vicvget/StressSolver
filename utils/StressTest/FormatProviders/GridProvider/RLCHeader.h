#ifndef RLC_HEADER_FILE_PROVIDER_H
#define RLC_HEADER_FILE_PROVIDER_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using std::string;
using std::ifstream;
using std::ofstream;

struct RLCHeader
{
	/*
	Версия формата RLC-файла.
	Всегда < 0, так как в первом формате сначала происходит
	запись количества точек сетки, а оно всегда положительное.
	Поэтому если версия больше нуля, то файл однозначно первого формата
	*/
	int _version;
	unsigned int _meshPos;
	unsigned int _normalsPos;
	unsigned int _freeSolverGridParamsPos;
	RLCHeader() : _version(0), _meshPos(0), _normalsPos(0), _freeSolverGridParamsPos(0) {};
};

class RLCHeaderFileProvider
{
public:
	RLCHeaderFileProvider() {};

public:
	static bool IsNewRLCFormat(const string &rlcFileName);
	static RLCHeader* LoadRLCHeader(ifstream &ifs);
	static RLCHeader* LoadRLCHeader(const string &rlcFileName);
	static void SaveRLCHeader(ofstream &ofs, RLCHeader* rlcHeader);
};

#endif