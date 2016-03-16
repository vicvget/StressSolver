#pragma once


#include "ParamsSetBase.h"
//#include "../../../../fcore/Calculator/Calculator.h"

#include <vector>
#include <string>


using std::vector;
using std::string;


#define SOLVER_TYPES_NUMBER 2
enum SolverTypes 
{
	ST_Thermal = 0,
	ST_StressStrain = 1
};


class SolverSpecialParamsTypeInfo
{
public:
	static size_t GetParamsCount(SolverTypes type);
	static void Init(SolverTypes type, vector<string> &params);
};


/** Специальные параметры решателя (например, свойства материала)
* класс для связывания наборов параметров для разных типов границ
*
* @author Getmanskiy Victor
*/
class SolverSpecialParams: public ParamsSetBase
{
	// тип граничных поерхностей
	SolverTypes _type;
	// номера граничных поверхностей
	vector<string> _params;
	// TODO: переделать на string и на инициализацию из String

public:

	SolverTypes GetType() const {return _type;}
	void SetType(SolverTypes val){_type = val;}
	//double* GetParamsPointer() {return _params.size() == 0 ? nullptr : &_params[0];}

	SolverSpecialParams
		(
			SolverTypes type,
			const string& name
		)
		:
			ParamsSetBase(name),
			_type(type)
	{
		Init();
	}

	SolverSpecialParams
		(
			SolverTypes type
		)
		:
			_type(type)
	{
		Init();
	}

	// только для wrapper-a!
	SolverSpecialParams()
		:
			_type()
	{
	}

	/** Заполняет набор параметров по умолчанию 
	* в зависимости от типа граничного условия
	*
	* @param SolverSpecialType - набор параметров по умолчанию 
	*/
	void CreateDefault(SolverTypes SolverType);
	
	size_t CountParams() const {return _params.size();};

	/** Установить значение параметра по индексу парметра
	* @param id - индекс параметра
	* @param newValue - новое значение параметра
	*/
	void SetParam(size_t id, string const& newValue);
	
	string GetParam(int id) const {return _params[id];}
	//double GetParam(Calc::Calculator* calc, int id) const {return calc->Eval(_params[id].c_str());}

#pragma region accessors
	vector<string>& GetParams(){return _params;}
	void SetParams(vector<string> val){_params = val;}
#pragma endregion

//TODO:
#pragma region overriden
	void Init() {SolverSpecialParamsTypeInfo::Init(_type,_params);};
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);
#pragma endregion

};

