#ifndef FFreeSolverH
#define FFreeSoverH

#include "../../../Fcore/fcore.h"
#include "MphParamsSet/FreeSolverGridParams.h"
#ifdef USE_PROVIDER_FRM_OLD
#	include "../../../Fcore/FrundFacade/MphParamsSet/SolverParamsBase.h"
#	include "../../../Fcore/FrundFacade/MphParamsSet/SolverSpecialParams.h"
#	include "../../../Fcore/FrundFacade/MphParamsSet/SolverIntParams.h"
#else
#include "MphParamsSet/SolverParamsBase.h"
#include "MphParamsSet/SolverSpecialParams.h"
#include "MphParamsSet/SolverIntParams.h"
#endif

#include <vector>
#include <string>
#include <iostream>
using std::vector;
using std::string;

///**
//*  ласс дл€ междисциплинарных свойств тела (набор дополнительных решателей)
//*
//* ¬ заголовке лежат все наборы параметров и выбранные значени€    
//*
//* @author Gromov Eugene
//*/
class FFreeSolver
{
private:
	// им€
	string _name;
	// индекс файла ( в нотации ‘–”Ќƒ 100000...)
	string _fileId;
	// номер
	int _number;
	// тип решател€
	SolverTypes _solverType;	
	// номер набора параметров интегрировани€ из заголовка
	size_t _intParamsId;
	// параметры сетки
	FreeSolverGridParams _gridParams;
	// номер набора специальных параметров решател€
	size_t _specialParamsId;
	// файл с сеткой - посмотреть, нужен ли он
	// используетс€ дл€ загрузки грид парамсов из rlc-файла
	// в дальнейшем он может помен€тьс€ и если каждый раз грузить из него, то можно нарватьс€ на изменение решател€ или на грид пармсы, относ€щиес€ к другой геометрии
	// свободный решатель полностью изолирован от тела, поэтому грид парамс€ должны быть в нем
	string _gridFileName;
	
	/* 
	   ѕор€док сохранени€ в файл:
	   1 - им€ решател€ _name
	   2 - тип решател€ _solverType
	   3 - им€ файла с сеткой _gridFileName
	   4 - активен ли решатель _isEnabled
	   5 - идентификатор параметра интегрировани€ _intParamsId
	   6 - идентификатор специального параметра _specialParamsId
	   7 - _readIco
	   8 - _writeIco
	   9 - грид парамсы
	*/
	
	bool _readIco;  // признак чтени€ ico
	bool _writeIco; // признак записи ico

	// признак активации решател€
	bool _isEnabled;

public:

	FFreeSolver();

	FFreeSolver(const char *fileName);
	FFreeSolver(string &name, string &fileName, SolverTypes solverType);

	void Init();

#pragma region Getters

	string GetFileIndex() const { return _fileId; }
	string FileId() const { return _fileId; }
	string Name() const { return _name; }
	int Number() const { return _number; }
	SolverTypes Type() const { return _solverType; }
	string GridFileName() const { return _gridFileName; }
	FreeSolverGridParams GridParams() const { return _gridParams; }
	bool IsEnabled() const { return _isEnabled; }
	int IntParamsId() const { return static_cast<int>(_intParamsId); }
	int SpecParamsId() const { return static_cast<int>(_specialParamsId); }
	bool ReadIco() const { return _readIco; }
	bool WriteIco() const { return _writeIco; }

#pragma endregion

#pragma region Setters

	void Name(string val) { _name = val; }
	void FileId(string val) { _fileId = val; }
	void Number(int val) { _number = val; }
	void Type(SolverTypes solverType) { _solverType = solverType; }
	void GridFileName(string gridFileName) { _gridFileName = gridFileName; }
	void GridParams(const FreeSolverGridParams &freeSolverGridParams) { _gridParams = freeSolverGridParams; }
	void IsEnabled(bool isEnabled) { _isEnabled = isEnabled; }
	void IntParamsId(int intParamsId) { _intParamsId = static_cast<size_t>(intParamsId); }
	void SpecParamsId(int specialParamsId) { _specialParamsId = static_cast<size_t>(specialParamsId); }
	void ReadIco(bool readIco) { _readIco = readIco; }
	void WriteIco(bool writeIco) { _writeIco = writeIco; }

#pragma endregion

	void Load(const char* path);
	int Load(ifstream& stream);
	void Save(ofstream& ofs) const;
	void SaveByIndex() const;
};

#endif