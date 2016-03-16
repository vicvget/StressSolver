#ifndef FSOLVE_PARAMS_SET_H

#define FSOLVE_PARAMS_SET_H


#include "FParamsSet.h"
#include "../../../fcore/fcore.h"
#include "BugStringsNoUTF.h"
#include "TransformMatrix.h"

#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <iomanip>


using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::pair;
using std::map;


typedef pair<int, TransformationMatrix> BodyNumberMatrixPair;
typedef map<int, TransformationMatrix> MatrixMap;


class IcoParams
{
private:

	MatrixMap _matrixMap;
public:
	string name;
	string icoFilePath;
	
	MatrixMap::iterator begin()
	{
		return _matrixMap.begin();
	}

	MatrixMap::const_iterator begin() const
	{
		return _matrixMap.begin();
	}

	MatrixMap::iterator end()
	{
		return _matrixMap.end();
	}

	MatrixMap::const_iterator end() const
	{
		return _matrixMap.end();
	}

	MatrixMap& GetMatrixMap()
	{
		return _matrixMap;
	}

	void SetMatrixMap
		(
			const MatrixMap& val
		)
	{
		_matrixMap = val;
	}

	void InsertIcoBodyMatrix
		(
			int bodyNumber,
			const TransformationMatrix& matrix
		)
	{
		_matrixMap[bodyNumber] = matrix;
	}

	bool GetIcoBodyMatrix
		(
			int bodyNumber,
			TransformationMatrix& matrix
		)	const
	{
		MatrixMap::const_iterator it;

		it = _matrixMap.find(bodyNumber);
		if (it == _matrixMap.end())
		{
			return false;
		}
		for (int i = 0; i < 3; i++)
		{
			matrix.ShiftVector[i] = it->second.ShiftVector[i];
		}
		for (int i = 0; i < 9; i++)
		{
			matrix.RotationMatrix[i] = it->second.RotationMatrix[i];
		}

		for (int i = 0; i < 9; i++)
		{
			matrix.SomeAngles[i] = it->second.SomeAngles[i];
		}

		return true;
	}

	bool GetIcoBodyMatrix
		(
			int bodyNumber,
			float (&matrix)[12]
		)	const
	{
		MatrixMap::const_iterator it;

		it = _matrixMap.find(bodyNumber);
		if (it == _matrixMap.end())
		{
			return false;
		}
		
		for (int i = 0; i < 3; i++)
		{
			matrix[i] = it->second.ShiftVector[i];
		}
		for (int i = 0; i < 9; i++)
		{
			matrix[i + 3] = it->second.RotationMatrix[i];
		}
		/*for (int i = 0; i < 9; i++)
		{
			matrix.SomeAngles[i] = it->second.SomeAngles[i];
		}*/

		return true;
	}

	IcoParams()
		:
			name(NAME_DEFAULT),
			icoFilePath(FILE_DEFAULT_INITIAL_CONDITIONS)
	{
		//name = new char[STRLIMIT];
		//icoFilePath = new char[STRLIMIT];
		//strcpy(name, NAME_DEFAULT);
		//strcpy(icoFilePath, FILE_DEFAULT_INITIAL_CONDITIONS);
	}
	~IcoParams()
	{
		/*if (name)
		delete name;

		if(icoFilePath)
		delete icoFilePath;*/
	}
	void Init()
	{
		name = NAME_DEFAULT;
		icoFilePath = FILE_DEFAULT_INITIAL_CONDITIONS;
		Load();
		//strcpy(name, NAME_DEFAULT);
		//strcpy(icoFilePath, FILE_DEFAULT_INITIAL_CONDITIONS);
	}

	void Load()
	{
		ifstream ifs(icoFilePath);

		if (ifs.is_open())
		{
			int bodyNumber;

			ifs	>> bodyNumber;
			while (!ifs.eof())
			{
				TransformationMatrix matrix;

				for (int i = 0; i < 3; i++)
				{
					ifs	>> matrix.ShiftVector[i];
				}


				for (int i = 0; i < 9; i++)
				{
					ifs	>> matrix.SomeAngles[i];
				}

				for (int i = 0; i < 9; i++)
				{
					ifs	>> matrix.RotationMatrix[i];
				}
				if (_matrixMap.find(bodyNumber) == _matrixMap.end())
				{
					_matrixMap.insert(BodyNumberMatrixPair(bodyNumber, matrix));
				}
				ifs	>> bodyNumber;
			}
			ifs.close();
		}
	}

	void Save() const
	{
		ofstream ofs(icoFilePath);

		ofs.setf(std::ios_base::uppercase | std::ios_base::scientific | std::ios_base::right);
		ofs.precision(6);
		if (ofs.is_open())
		{
			MatrixMap::const_iterator it;

			for (it = _matrixMap.begin(); it != _matrixMap.end(); ++it)
			{
				int bodyNumber  = (*it).first;
				TransformationMatrix matrix;
				if(!GetIcoBodyMatrix(bodyNumber,matrix))
					continue;
				ofs << "           " << bodyNumber << std::endl;

				for (int i = 0; i < 3; i++)
				{
					ofs << std::setw(16) << matrix.ShiftVector[i] << ' ';
				}
				for (int i = 0; i < 9; i++)
				{
					ofs << std::setw(16) << matrix.SomeAngles[i];
					if(i == 0 || i == 2 || i == 6 || i == 8)
						ofs	<< std::endl;
					else
						ofs	<< ' ';
				}
				for (int i = 0; i < 9; i++)
				{
					ofs << std::setw(16) << matrix.RotationMatrix[i];
					if(i == 3 || i == 7 || i == 8)
						ofs	<< std::endl;
					else
						ofs	<< ' ';
				}
			}//it
			ofs.close();
		}//ofs.is_open
	}//Save();
};

class ControlParams
{
public:

	string name;
	size_t numberOfSolveVariants;
	double velocity;
	double velocityStep;
	double solverStep;
	double solverOutputStep;
	size_t solveType;
	double time;
	double outputFromTime;
	size_t outputCode;

	ControlParams()
	{
		Init();
	}
	~ControlParams()
	{
	}


	void Init(bool isAssembly = false)
	{
		if(isAssembly)
		{
			name = STRING_FIX;
			numberOfSolveVariants=1;
			velocity=1;
			velocityStep=1;
			solverStep=0.00001;
			solverOutputStep=0.001;
			solveType=8;
			time=5;
			outputFromTime=0;
			outputCode=1;
		}
		else
		{
			name = NAME_DEFAULT;
			numberOfSolveVariants=1;
			velocity=1;
			velocityStep=1;
			solverStep=0.00001;
			solverOutputStep=0.001;
			solveType=3;
			time=5;
			outputFromTime=0;
			outputCode=1;
		}
	}
};

class ExtraControlParams
{
public:
	string name;
	double factorOfNondissipativeRegulateForces;
	double factorOfDissipativeRegulateForces;
	double maxResidualInJoints;
	double factorOfAMTFrequency; //AdjacencyMatrixTriangulationFrequency
	double factorOfResidualRegulation;
	double constraintsApplyStep;

	ExtraControlParams()
	{
		Init();
		//name = NAME_DEFAULT;
		//name = new char[STRLIMIT];
		//strcpy(name, NAME_DEFAULT);
	}
	~ExtraControlParams()
	{

		/*if(name)
		delete name;*/
	}

	void Init(bool isAssembly = false)
	{
		if(isAssembly)
		{
			name = STRING_FIX;
			//strcpy(name, NAME_DEFAULT);
			factorOfNondissipativeRegulateForces = 500000000;
			factorOfDissipativeRegulateForces	 = 2000;
			maxResidualInJoints					 = 0.000001;
			factorOfAMTFrequency				 = 1; //AdjacencyMatrixTriangulationFrequency
			factorOfResidualRegulation			 = 250;
			constraintsApplyStep				 = 1;
		}
		else
		{
			name = NAME_DEFAULT;
			//strcpy(name, NAME_DEFAULT);
			factorOfNondissipativeRegulateForces = 0;
			factorOfDissipativeRegulateForces	 = 0;
			maxResidualInJoints					 = 10.001;
			factorOfAMTFrequency				 = 1; //AdjacencyMatrixTriangulationFrequency
			factorOfResidualRegulation			 = 250;
			constraintsApplyStep				 = 1;
		}
	}
};

class SpecialParams
{
public:
	string name;
	string leftTrackFilePath;
	string rightTrackFilePath;
	size_t trackType; // 0 = none, 1 = single, 2 = sin, 3 = random
	bool useInitialVelocity;
	bool direction; // true = 'XPlus', false = 'XMinus'
	double leftTrackParams[6];
	double rightTrackParams[6];
	bool assumeAllJointsFlexible;

	SpecialParams()
	{
		Init();
		/*name = new char[STRLIMIT];
		leftTrackFilePath = new char[STRLIMIT];
		rightTrackFilePath = new char[STRLIMIT]; 

		const char* FILE_TRACK_PATH_DEFAULT = "none";
		strcpy(name, NAME_DEFAULT);
		strcpy(leftTrackFilePath, FILE_TRACK_PATH_DEFAULT);
		strcpy(rightTrackFilePath, FILE_TRACK_PATH_DEFAULT);*/

	}
	~SpecialParams()
	{
		/*if (name)
		delete name;
		if(leftTrackFilePath)
		delete leftTrackFilePath;
		if (rightTrackParams)
		delete rightTrackParams;*/

	}
	void Init()
	{
		const char* FILE_TRACK_PATH_DEFAULT = "none";
		name = NAME_DEFAULT;
		//strcpy(name, NAME_DEFAULT);
		trackType = 1; //starting from 1 = NO
		useInitialVelocity = 0;
		direction = true; 
		leftTrackFilePath = FILE_TRACK_PATH_DEFAULT;
		rightTrackFilePath = FILE_TRACK_PATH_DEFAULT;
		//strcpy(leftTrackFilePath, FILE_TRACK_PATH_DEFAULT);
		//strcpy(rightTrackFilePath, FILE_TRACK_PATH_DEFAULT);
		leftTrackParams[0] = rightTrackParams[0] = 0;
		leftTrackParams[1] = rightTrackParams[1] = 0.5;
		leftTrackParams[2] = rightTrackParams[2] = 0.5;
		leftTrackParams[3] = rightTrackParams[3] = 0.5;
		leftTrackParams[4] = rightTrackParams[4] = 1.0;
		leftTrackParams[5] = rightTrackParams[5] = 0;
		assumeAllJointsFlexible = false;
	}
};

class ParamSetIds
{
public:
	string name;
	size_t IcoParamsId;
	size_t ControlParamsId;
	size_t ExtraControlParamsId;
	size_t SpecialParamsId;
	size_t VariableParamsId;


	ParamSetIds()
	{
		Init();
	}
	~ParamSetIds()
	{

	}
	void Init(bool isAssembly = false)
	{
		if(isAssembly)
		{
			name = STRING_FIX;
			IcoParamsId = 0;
			ControlParamsId = 1;
			ExtraControlParamsId = 1;
			SpecialParamsId = 0;
			VariableParamsId = 0;
		}		
		else
		{
			name = NAME_DEFAULT;
			IcoParamsId = 0;
			ControlParamsId = 0;
			ExtraControlParamsId = 0;
			SpecialParamsId = 0;
			VariableParamsId = 0;
		}
	}
};

class FSolveParamsSet: public FParamsSet
{
	vector<IcoParams>			icoParamsSets;
	vector<ControlParams>		controlParamsSets;
	vector<ExtraControlParams>	extraControlParamsSets;
	vector<SpecialParams>		specialParamsSets;
	vector<ParamSetIds>			paramSetIdsSets;
	void OutputIcoParams(ofstream& ofs, size_t id) const;
public:
	FSolveParamsSet();
	~FSolveParamsSet();

	void AddIcoParamsSet(const IcoParams& params) {icoParamsSets.push_back(params);};
	void AddControlParamsSet(const ControlParams& params) {controlParamsSets.push_back(params);};
	void AddExtraControlParamsSet(const ExtraControlParams& params) {extraControlParamsSets.push_back(params);};
	void AddSpecialParamsSet(const SpecialParams& params) {specialParamsSets.push_back(params);};
	void AddParamSetIdsSet(const ParamSetIds& params) {paramSetIdsSets.push_back(params);};
	void FullClean();

	size_t CurrentId() const { return currentId; }
	void CurrentId(size_t val) { currentId = val; }

	size_t Size() const {return paramSetIdsSets.size();}

	//public bool Select(int id) {currentSetId = id;}

	/** Считывание параметров из входного потока
	* @param ifs - входной поток
	*/
	void Read(ifstream& ifs);

	/** Запись параметров в выходной поток
	* @param ofs - выходной поток
	*/
	void Write(ofstream& ofs) const;

	/** Вывод в uprf Параметров расчета
	* @param ofs - выходной поток
	*/
	void OutputControlParams(ofstream& ofs) const;

	/** Вывод в way.cnt имен файлов с левым и правым профилем дороги
	* @param ofs - выходной поток
	*/
	void OutputRoadParams(ofstream& ofs) const;

	/** Вывод в epsilon.dat дополнительных отладочных параметров
	* @param ofs - выходной поток
	*/
	void OutputExtraParams(ofstream& ofs) const;

	/** Вывод в default.ico начальных условий
	* @param ofs - выходной поток
	*/
	void OutputIcoParams(ofstream& ofs) const;

	/** Вывод ico по умолчанию
	*/
	void OutputDefaultIcoParams() const;
	// Доступ к наборам переменных
#pragma region Accessors

	ParamSetIds* GetParamSetIdsSet(int id = -1)
	{
		return id == -1 ? &paramSetIdsSets.at(currentId) : &paramSetIdsSets.at(id);
	}
	
	const ParamSetIds* GetParamSetIdsSet(int id = -1) const
	{
		return id == -1 ? &paramSetIdsSets.at(currentId) : &paramSetIdsSets.at(id);
	}

	IcoParams*			GetIcoParamsSet(int id = -1)		  {return id == -1 ? &icoParamsSets.at(GetParamSetIdsSet()->IcoParamsId) : &icoParamsSets.at(id);}
	ControlParams*		GetControlParamsSet(int id = -1)	  {return id == -1 ? &controlParamsSets.at(GetParamSetIdsSet()->ControlParamsId) : &controlParamsSets.at(id);}
	ExtraControlParams* GetExtraControlParamsSet(int id = -1) {return id == -1 ? &extraControlParamsSets.at(GetParamSetIdsSet()->ExtraControlParamsId): &extraControlParamsSets.at(id);}

	const IcoParams&			GetIcoParamsSet(int id = -1)			const { return id == -1 ? icoParamsSets.at(GetParamSetIdsSet()->IcoParamsId) : icoParamsSets.at(id); }
	const ControlParams&		GetControlParamsSet(int id = -1)		const { return id == -1 ? controlParamsSets.at(GetParamSetIdsSet()->ControlParamsId) : controlParamsSets.at(id); }
	const ExtraControlParams&   GetExtraControlParamsSet(int id = -1)	const { return id == -1 ? extraControlParamsSets.at(GetParamSetIdsSet()->ExtraControlParamsId) : extraControlParamsSets.at(id); }
	const ParamSetIds&			GetParamSetIdsSetRef(int id = -1)		const { return id == -1 ? paramSetIdsSets.at(currentId) : paramSetIdsSets.at(id); }

	SpecialParams* GetSpecialParamsSet(int id = -1)
	{
		return id == -1 ? &specialParamsSets.at(GetParamSetIdsSet()->SpecialParamsId) : &specialParamsSets.at(id);
	}

	const SpecialParams* GetSpecialParamsSet(int id = -1) const
	{
		return id == -1 ? &specialParamsSets.at(GetParamSetIdsSet()->SpecialParamsId) : &specialParamsSets.at(id);
	}

	vector<ParamSetIds> ParamSetIdsSets() const { return paramSetIdsSets; }
	void ParamSetIdsSets(vector<ParamSetIds> val) { paramSetIdsSets = val; }

	vector<ControlParams> ControlParamsSets() const { return controlParamsSets; }
	void ControlParamsSets(vector<ControlParams> val) { controlParamsSets = val; }

	vector<IcoParams> IcoParamsSets() const { return icoParamsSets; }
	void IcoParamsSets(vector<IcoParams> val) { icoParamsSets = val; }

	vector<SpecialParams> SpecialParamsSets() const { return specialParamsSets; }
	void SpecialParamsSets(vector<SpecialParams> val) { specialParamsSets = val; }

	vector<ExtraControlParams> ExtraControlParamsSets() const { return extraControlParamsSets; }
	void ExtraControlParamsSets(vector<ExtraControlParams> val) { extraControlParamsSets = val; }


#pragma endregion

#pragma region Marshaling
	size_t GetIcoParamsCount()			const {return icoParamsSets.size();}
	size_t GetControlParamsCount()		const {return controlParamsSets.size();}
	size_t GetExtraControlParamsCount() const {return extraControlParamsSets.size();}
	size_t GetSpecialParamsCount()		const {return specialParamsSets.size();}
	size_t GetParamSetIdsSetsCount()	const {return paramSetIdsSets.size();}

	IcoParams*			GetIcoParamsSets()			{return &icoParamsSets[0];}
	ControlParams*		GetControlParamsSets()		{return &controlParamsSets[0];}
	ExtraControlParams* GetExtraControlParamsSets() {return &extraControlParamsSets[0];}
	SpecialParams*		GetSpecialParamsSets()		{return &specialParamsSets[0];}
	ParamSetIds*		GetParamSetIdsSets()		{return &paramSetIdsSets[0];}
#pragma endregion

};


#endif // FSOLVE_PARAMS_SET_H