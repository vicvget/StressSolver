#include "SolverSpecialParams.h"
#include "../../../../fcore/wrappers/FileRoutines.h"

void SolverSpecialParams::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	fs::WriteCommentedLine(ofs, GetType(), "special type");
	fs::WriteCommentedLine(ofs, CountParams(), "count params");
	for(int i = 0; i < CountParams(); i++ )
	{
		ofs << _params[i] << std::endl;
	}
}

void SolverSpecialParams::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	size_t tmpCount;
	int tmpType;
	fs::ReadLineValue(ifs, tmpType);
	_type = (SolverTypes)tmpType;
	Init();
	fs::ReadLineValue(ifs, tmpCount);
	if(tmpCount < SolverSpecialParamsTypeInfo::GetParamsCount(_type))
	{
		// это на случай, если при изменении формата надо прочиать файл со старыми даными
		_params.resize(SolverSpecialParamsTypeInfo::GetParamsCount(_type));
	}
	else
		_params.resize(tmpCount);
	for(int i = 0; i < tmpCount; i++)
	{
		fs::ReadLineValue(ifs, _params[i]);
	}
}

/** Заполняет набор параметров по умолчанию 
* в зависимости от типа граничного условия
*
* @param SolverType - набор параметров по умолчанию 
*/
void SolverSpecialParams::CreateDefault(SolverTypes solverType)
{
	SolverSpecialParamsTypeInfo::Init(solverType, _params);
}
void SolverSpecialParamsTypeInfo::Init(SolverTypes type, vector<string> &params)
{
	params.resize(GetParamsCount(type));
	switch (type)
	{
	case ST_Thermal:
		/**
		* [0] TPR - теплопроводность;
		* [1] TPL - теплоемкость;
		* [2] PLOT - плотность;				
		* [3] T0 - начальная температура
		* [4] коэффициент масштабирования теплового потока (умножается на tfunc);
		*/
		params[0] = "47";
		params[1] = "462";
		params[2] = "7800";
		params[3] = "20";
		params[4] = "1";
		break;

	case ST_StressStrain:
		/**
		* E - модуль упругости;
		* DM - коэффициент демпфирования
		* Pho - плотность
		* Sf - коэффициент масштабирования
		* S - флаг статического расчета (0 - динамический расчет, 1 - статический)
		*/
		params[0] = "1e5";
		params[1] = "50";
		params[2] = "7800";
		params[3] = "1";
		params[4] = "0";
		break;
	}
}

size_t SolverSpecialParamsTypeInfo::GetParamsCount(SolverTypes type)
{
	switch (type)
	{
	case ST_Thermal:
		return 5;

	case ST_StressStrain:
		return 5;
	}
	return 0;
}
void SolverSpecialParams::SetParam(size_t id, string const&newValue)
{
	if (id < _params.size())
		_params[id] = newValue;
}

