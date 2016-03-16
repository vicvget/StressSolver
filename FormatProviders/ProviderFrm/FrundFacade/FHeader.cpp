#include "FHeader.h"

#include "../../../fcore/wrappers/FileRoutines.h"

#include <fstream>
#include <string>


using std::endl;


FHeader::FHeader()
	:
		_taskName("newtask"),
		_taskCode("newmodel"),
		_mphHeader(nullptr)
{
	_fileId = "30000001";
//	this->SetSolveParamsSet(0);
//	this->SetVariablesParamsSet(0);
}

// Загрузка заголовка модели ФРУНД из файла
bool FHeader::Load(ifstream& stream)
{
	fs::ReadLineString(stream, _taskName);
	fs::ReadLineString(stream, _taskCode);
	ReadVariables(stream);
	ReadMacros(stream);
	ReadSolverParams(stream);
	ReadMeasureParams(stream);
	//if(stream.eof()) return 0;

	return true;
}

void FHeader::Load(const char* path)
{
	DefaultLoad(path);
	
	string strPath(path);
	int pos = strPath.find_first_of(EXT_MODEL_ELEMENT);
	
	if(pos == string::npos)
		exceptions::ThrowFileInvalidExtension(strPath);

	_fileId = strPath.substr(0, pos);
	_mphHeader = nullptr;

	string fileMph = _fileId + EXT_MULTIPHYSICS;
	ifstream streamGeo(fileMph);
	if(streamGeo.is_open())
	{
		_mphHeader = new FMphHeader(fileMph.c_str());
		streamGeo.close();
	}
}

void FHeader::Save(ofstream& ofs) const
{
	ofs << _taskName << endl 
		<< _taskCode << endl;
	WriteVariables(ofs);
	WriteMacros(ofs);
	WriteSolverParams(ofs);
	WriteMeasureParams(ofs);
	// TODO: разобраться, что это за данные
	ofs << "0 0 0 0 2 0.000000 100.000000 0.010000\n0.000000 0 Reserved string\n";
}

void FHeader::SaveByIndex() const
{
	FElement::SaveByIndex();
	
	if(_mphHeader != nullptr)
	{
		string path = _fileId + EXT_MULTIPHYSICS;
		ofstream ofs(path);
		if(ofs.is_open())
		{
			_mphHeader->Save(ofs);
			ofs.close();
		}
		else
		{
			exceptions::ThrowFileNotOpened(path);
		}
	}
}


// Загрузка переменных из файлового потока
void FHeader::ReadVariables(ifstream &stream)
{
	variablesGroupsSet.Read(stream);
}


void FHeader::WriteVariables(ofstream &ofs) const
{
	variablesGroupsSet.Write(ofs);
}

void FHeader::WriteMacros(ofstream &ofs) const
{
	int strCount = 0;
	for(int i = 0; i < _macros.size(); i++)
		strCount += _macros[i].StrCount();
	ofs << strCount << endl;
	for(int i = 0; i < _macros.size(); i++)
	{
		_macros[i].Write(ofs);
	}	
}

void FHeader::WriteSolverParams(ofstream &ofs) const
{
	solveParamsSet.Write(ofs);
}

void FHeader::WriteMeasureParams(ofstream &ofs) const
{
	this->measureParamsSet.Write(ofs);
}

void FHeader::ReadMacros(ifstream &stream)
{
	int strCount;
	fs::ReadLineValue(stream,strCount);
	while(strCount > 0)
	{
		FMacro tmpMacro(stream);
		strCount -= tmpMacro.StrCount();
		_macros.push_back(tmpMacro);
	}
}

void FHeader::ReadSolverParams(ifstream &stream)
{
	solveParamsSet.Read(stream);
}

void FHeader::ReadMeasureParams(ifstream &stream)
{
	measureParamsSet.Read(stream);
}


// Вывод в текстовый файловый поток переменных текущего набора
void FHeader::OutputVariables(ofstream &stream, int isRType) const
{
	if(isRType)
		variablesGroupsSet.CurrentVariablesGroup().Output(stream);
}

// Вывод в текстовый файловый поток макросов
void FHeader::OutputMacros(ofstream &stream, int isRType) const
{
	for(int i = 0; i < _macros.size(); i++)
	{
		_macros[i].Output(stream, isRType);
	}
}

// Вывод в текстовый файловый поток заголовка расчетной схемы
void FHeader::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	const char* stpr = "STPR ";
	const char* stpm = "STPM ";
	const char* prefix = isRType ? stpr : stpm;
	//char* prefix = "STPM ";

	//if(isRType) prefix[3] = 'R';
	stream << "LIN1 " << prefix << '\'' << _taskName << "\'\n";
	stream << "LIN2 " << '\'' << _taskCode << "\'\n";
}


//// Вывод управляющих параметров для расчета(uprf,epsilon.dat,default.ico)
//void FHeader::OutputSolveControlParams()
//{
//	solveParamsSet.Output(stream);
//}

// Вывод управляющих параметров для расчета(UPRF)
void FHeader::OutputSolveControlParams(ofstream &stream) const
{
	solveParamsSet.OutputControlParams(stream);
}

// Вывод управляющих параметров для анализа (UPRF)
void FHeader::OutputAnazControlParams(ofstream &stream) const
{
	measureParamsSet.Output(stream);
}

/** Формирует управляющие файлы для шага расчет
* @param isRoadExist - признак наличия дороги в модели
*/
void FHeader::OutputSolveParams(bool isRoadExist) const
{
	ofstream ofs;
	
	if (IsTireExistInMacros())
	{
		ofs.open(FILE_TIRES_PATHS);
		if(ofs.is_open())
		{
			// TODO: не работает, отладить!
			//ofs << fs::CombinePathEnv(string(DIR_DATA),string(FILE_TIRE_DEFAULT_PARAMS_1)) << std::endl
			//	<< fs::CombinePathEnv(string(DIR_DATA),string(FILE_TIRE_DEFAULT_PARAMS_2)) << std::endl;
			// копирование файлов "tfy2123.dat" и "tmz2123.dat" из директории "Data" в каталог с моделью
			if (!fs::FileExist(FILE_TIRE_DEFAULT_PARAMS_1))
			{
				fs::CopyFileToFolder
					(
						fs::CombinePathEnv
							(
								DIR_DATA,
								FILE_TIRE_DEFAULT_PARAMS_1
							),
						fs::GetCurrentDir()
					);
			}
			if (!fs::FileExist(FILE_TIRE_DEFAULT_PARAMS_2))
			{
				fs::CopyFileToFolder
					(
						fs::CombinePathEnv
							(
								DIR_DATA,
								FILE_TIRE_DEFAULT_PARAMS_2
							),
						fs::GetCurrentDir()
					);
			}
			ofs << string(FILE_TIRE_DEFAULT_PARAMS_1) << std::endl
				<< string(FILE_TIRE_DEFAULT_PARAMS_2) << std::endl;
			ofs.close();
			ofs.clear();
		}
		else
			exceptions::ThrowFileNotOpened(FILE_TIRES_PATHS);
	}
	if(isRoadExist)
	{
		ofs.open(FILE_ROAD_PROFILES_PATHS);
		if(ofs.is_open())
		{
			solveParamsSet.OutputRoadParams(ofs);
			ofs.close();
			ofs.clear();
		}
		else
			exceptions::ThrowFileNotOpened(FILE_ROAD_PROFILES_PATHS);
	}

	ofs.open(FILE_CONTROL_PARAMS);
	if(ofs.is_open())
	{
		solveParamsSet.OutputControlParams(ofs);
		ofs.close();
		ofs.clear();
	}
	else
		exceptions::ThrowFileNotOpened(FILE_CONTROL_PARAMS);
	ofs.open(FILE_EXTRA_CONTROL_PARAMS);
	if(ofs.is_open())
	{
		solveParamsSet.OutputExtraParams(ofs);
		ofs.close();
		ofs.clear();
	}
	else
		exceptions::ThrowFileNotOpened(FILE_EXTRA_CONTROL_PARAMS);
	
	ofs.open(FILE_INITIAL_CONDITIONS);
	if(ofs.is_open())
	{
		solveParamsSet.OutputIcoParams(ofs);
		//solveParamsSet.OutputDefaultIcoParams()
		ofs.close();
		ofs.clear();
	}
	else
		exceptions::ThrowFileNotOpened(FILE_INITIAL_CONDITIONS);
}

/** Выбрать набор управляющих параметров для расчета
* @param paramsSetId - номер набора параметров
*/
void FHeader::SetSolveParamsSet(int paramsSetId)
{
	solveParamsSet.Select(paramsSetId);
}

/** Выбрать набор управляющих параметров для построения графиков
* @param paramsSetId - номер набора параметров
*/
void FHeader::SetPostprocessorControlParamsSet(int paramsSetId)
{
	measureParamsSet.Select(paramsSetId);
}
	
/** Выбрать набор параметров модели
* @param paramsSetId - номер набора параметров
*/
void FHeader::SetVariablesParamsSet(int paramsSetId)
{
	// если -1, то оставить текущий загруженный по умолчанию
	if(paramsSetId != -1)
		variablesGroupsSet.Select(paramsSetId);
}

FHeader::~FHeader()
{
	if (_mphHeader != nullptr)
		delete _mphHeader;
}

bool FHeader::IsTireExistInMacros() const
{
	for (size_t i = 0; i < _macros.size(); i++)
	{
		if (_macros[i].Type().find("tire") != string::npos)
		{
			return true;
		}
	}

	return false;
}
