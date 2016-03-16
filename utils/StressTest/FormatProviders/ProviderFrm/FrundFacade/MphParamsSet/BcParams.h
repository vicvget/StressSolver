#ifndef BC_PARAMS_H

#define BC_PARAMS_H


#include "ParamsSetBase.h"
#include "../../../../fcore/Calculator/Calculator.h"

#include <vector>
#include <string>


using std::vector;
using std::string;


#define BC_TYPES_NUMBER 5

enum BcTypes
{
	BCT_ThermalBoundary1 = 0,
	BCT_ThermalBoundary2 = 1,
	BCT_ThermalBoundary3 = 2,
	BCT_StressStrainBoundaryForce = 3,
	BCT_StressStrainBoundarySealing = 4,
};

class BcTypeInfo
{
public:
	static size_t GetParamsCount(BcTypes type);
	static int GetCouplerId(BcTypes type);
	static void Init(BcTypes type, vector<string> &params);
};

/** Параметры граничных условий
* класс для связывания наборов параметров для разных типов границ
* с объектами BcSurface
*
* @author Getmanskiy Victor
*/
class BcParams: public ParamsSetBase
{
	// тип граничных поерхностей
	BcTypes _type;
	// номера граничных поверхностей
	vector<string> _params;

public:

	//int GetType(){return (int)_type;}
	//double* GetParamsPointer() {return _params.size() == 0 ? nullptr: &_params[0];}

	BcParams
		(
			BcTypes type,
			const string& name
		)
		:
			ParamsSetBase(name),
			_type(type)
	{
		Init();
	}

	BcParams
		(
			BcTypes type
		)
		:
			_type(type)
	{
		Init();
	}

	 // только для wrapper-a!
	BcParams()
		:
			_type()
	{
	}

	/** Заполняет набор параметров по умолчанию 
	* в зависимости от типа граничного условия
	*
	* @param bcType - набор параметров по умолчанию 
	*/
	void CreateDefault(BcTypes bcType);
	
	size_t CountParams() const {return _params.size();};

#pragma region accessors
	vector<string>& Params(){return _params;}
	vector<string> GetParams() const {return _params;}
	void SetParams(vector<string> val) { _params = val; }

	BcTypes GetType() const {return _type;}
	void SetType(BcTypes val) {_type = val;}
#pragma endregion

//TODO:
#pragma region overriden	
	void Init() {BcTypeInfo::Init(_type,_params);};
	void Save(ofstream& ofs) const;
	void Load(ifstream& ifs);
#pragma endregion

};


#endif // BC_PARAMS_H