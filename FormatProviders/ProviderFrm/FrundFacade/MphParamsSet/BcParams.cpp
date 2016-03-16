#include "BcParams.h"

#include "../../../../fcore/wrappers/FileRoutines.h"

#include <algorithm>


using std::find;


void BcParams::Save(ofstream& ofs) const
{
	DefaultSave(ofs);
	//ofs << GetType() << std::endl;
	//ofs << CountParams() << std::endl; // TODO: это убрать ?
	fs::WriteCommentedLine(ofs, GetType(), "bc type");
	fs::WriteCommentedLine(ofs, CountParams(), "count params");
	for(int i = 0; i < CountParams(); i++ )
	{
		ofs << _params[i] << std::endl;
	}
}

void BcParams::Load(ifstream& ifs)
{
	DefaultLoad(ifs);
	size_t tmpCount;
	int tmpType;
	fs::ReadLineValue(ifs, tmpType);
	_type = (BcTypes)tmpType;
	Init();
	fs::ReadLineValue(ifs, tmpCount);
	// tmpCount = BcTypeInfo::GetParamsCount(_type);

	if(tmpCount < BcTypeInfo::GetParamsCount(_type))
	{
		// это на случай, если при изменении формата надо прочиать файл со старыми даными
		// такого случая не будет!
		_params.resize(BcTypeInfo::GetParamsCount(_type));
		for(int i = tmpCount-1; i < BcTypeInfo::GetParamsCount(_type); i++)
			_params[i]="0";
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
* @param bcType - набор параметров по умолчанию 
*/
void BcParams::CreateDefault(BcTypes bcType)
{
	BcTypeInfo::Init(bcType, _params);
}
void BcTypeInfo::Init(BcTypes type, vector<string> &params)
{
	params.resize(GetParamsCount(type));
	switch (type)
	{
	case BCT_ThermalBoundary1:
		params[0] = "100";
		break;
	case BCT_ThermalBoundary2:
		params[0] = "100";
		params[1] = "0";
		break;
	case BCT_ThermalBoundary3:
		params[0] = "100";
		params[1] = "30";
		break;
	case BCT_StressStrainBoundaryForce:
		params[0] = "0";
		params[1] = "0";
		params[2] = "0";
		params[3] = "0";
		params[4] = "0";
		params[5] = "0";
		params[6] = "0";
		break;
	case BCT_StressStrainBoundarySealing:
		params[0] = "-1"; // '-' = sealed, '+' = unsealed
		params[1] = "-1";
		params[2] = "-1";
		params[3] = "-1";
		params[4] = "-1";
		params[5] = "-1";
		break;
	}
}

size_t BcTypeInfo::GetParamsCount(BcTypes type)
{
	switch (type)
	{
		case BCT_ThermalBoundary1:
			return 1;
		case BCT_ThermalBoundary2:
			return 2;
		case BCT_ThermalBoundary3:
			return 2;
		case BCT_StressStrainBoundaryForce:
			return 7;
		case BCT_StressStrainBoundarySealing:
			return 6;
		default: 
			return 0;
	}
	return 0;
}

int BcTypeInfo::GetCouplerId(BcTypes type)
{
	int id = -1;
	switch (type)
	{
		case BCT_ThermalBoundary2:
			id = 1; // второй параметр - номер СЕ
			break;
		case BCT_StressStrainBoundaryForce:
			id = 6; // шестой параметр - номер СЕ
			break;
		default:
			id = -1;
	}
	return id;
}