#include "FMeasureParamsSet.h"

#include "../../../fcore/wrappers/FileRoutines.h"
#include "../../../fcore/Exceptions/fcExceptions.h"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>
#include <functional>


using std::pair;
using std::binary_function;
using std::setw;
using std::endl;


// struct MeasureParams

MeasureParams::MeasureParams()
	:
		satellite(nullptr),
		name("default"),
		object(),
		subobject(),
		outputClass(),
		outputCode(),
		dof(),
		characteristic(),
		isLogarifmic(),
		params(),
		object1(),
		subobject1(),
		windowType(),
		fftBlockLength(),
		isMovingCS(),
		measureGroupId()
{
}

MeasureParams::MeasureParams
	(
		const MeasureParams* mp
	)
	:
		satellite(nullptr),
		object(),
		subobject(),
		outputClass(),
		outputCode(),
		dof(),
		characteristic(),
		isLogarifmic(),
		params(),
		object1(),
		subobject1(),
		windowType(),
		fftBlockLength(),
		isMovingCS(),
		measureGroupId()
{
	Copy(*mp);
}

MeasureParams::MeasureParams
	(
		const MeasureParams& mp
	)
	:
		satellite(nullptr),
		object(),
		subobject(),
		outputClass(),
		outputCode(),
		dof(),
		characteristic(),
		isLogarifmic(),
		params(),
		object1(),
		subobject1(),
		windowType(),
		fftBlockLength(),
		isMovingCS(),
		measureGroupId()
{
	Copy(mp);
}

MeasureParams::~MeasureParams()
{
	delete satellite;
}

MeasureParams& MeasureParams::operator = (const MeasureParams& mp)
{
#ifdef _DEBUG
	printf("MeasureParams& operator = (const MeasureParams& mp)\n");
#endif
	if (&mp == this)
	{
		// TODO: вставить выброс исключения

		return *this;
	}

	if (satellite != nullptr)
	{
		delete satellite;
		satellite = nullptr;
	}
	Copy(mp);

	return *this;
}

void MeasureParams::Copy(const MeasureParams& mp)
{
	CopyMembers(mp);
	if(mp.satellite != nullptr)
	{
	#ifdef _DEBUG
		printf("mp.satellite != nullptr\n");
	#endif
		satellite = new MeasureParams(mp.satellite);
	}
	else
	{
		satellite = nullptr;
	}
}

void MeasureParams::CopyMembers(const MeasureParams& mp)
{
	//memcpy(this, &mp, sizeof(MeasureParams));
	satellite = mp.satellite;
	name = mp.name;
	object = mp.object;
	subobject = mp.subobject;
	outputClass = mp.outputClass;
	outputCode = mp.outputCode;
	for (int i = 0; i < 6; i++)
	{
		dof[i] = mp.dof[i];
	}
	for (int i = 0; i < 3; i++)
	{
		characteristic[i] = mp.characteristic[i];
	}
	isLogarifmic = mp.isLogarifmic;
	for (int i = 0; i < 3; i++)
	{
		params[i] = mp.params[i];
	}
	object1 = mp.object1;
	subobject1 = mp.subobject1;
	windowType = mp.windowType;
	fftBlockLength = mp.fftBlockLength;
	isMovingCS = mp.isMovingCS;
	measureGroupId = mp.measureGroupId;
}


FMeasureParamsSet::FMeasureParamsSet()
{

}

FMeasureParamsSet::~FMeasureParamsSet()
{

}

bool ReadMeasure(ifstream& ifs, MeasureParams* it, string& groupName)
{
	stringstream str;
	bool res = true;

	res &= fs::ReadLineStringNoException(ifs, it->name);
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->object >> it->subobject >> it->outputClass >> it->outputCode;
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->dof[0] >> it->dof[1] >> it->dof[2] >> it->dof[3] >> it->dof[4] 
		>> it->dof[5];
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->characteristic[0] >> it->characteristic[1] >> it->characteristic[2];
	res &= fs::ReadLineValueNoException(ifs, it->isLogarifmic);
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->params[0] >> it->params[1] >> it->params[2];
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->object1 >> it->subobject1;
	res &= fs::ReadLineNoException(ifs, str);
	str >> it->windowType >> it->fftBlockLength >> it->isMovingCS;
	res &= fs::ReadLineStringNoException(ifs, groupName);

	return res;
}

void WriteMeasure(ofstream& ofs, const MeasureParams* it, const string& groupName)
{
	ofs << it->name << std::endl;
	ofs << it->object << ' ' << it->subobject << ' ' 
		<< it->outputClass << ' ' << it->outputCode << std::endl;
	for(int i = 0; i < 6; i++)
		ofs << it->dof[i] << ((i == 5) ? '\n' : ' ');
	for(int i = 0; i < 3; i++)
		ofs << it->characteristic[i] << ((i == 2) ? '\n' : ' ');
	ofs << it->isLogarifmic << std::endl;
	for(int i = 0; i < 3; i++)
		ofs << it->params[i] << ((i == 2) ? '\n' : ' ');
	ofs << it->object1 << ' '<< it->subobject1 << std::endl;
	ofs << it->windowType << ' ' << it->fftBlockLength << ' ' << it->isMovingCS << std::endl;
	ofs << groupName << std::endl;
}

void FMeasureParamsSet::Read(ifstream& ifs)
{
	size_t count;
	fs::ReadLineValue(ifs, count);
	currentId = 1;
	if(count != 0)
	{
		measureParamsSets.resize(count);
		vector<MeasureParams>::iterator it = measureParamsSets.begin();
		int measureParamsSetsId = 0;
		int measureGroupId = 0;
		int counter = 0;
		while(it != measureParamsSets.end())
		{
			string groupName, groupNameSatellite;
			if(!ReadMeasure(ifs, &(*it), groupName))
			{
				measureParamsSets.resize(counter);
				return;
			}
			counter++;
			if(it->object1 != 0)
			{
				// TODO: сделать генерацию UPRF для VIV satellite
				it->satellite = new MeasureParams();
				ReadMeasure(ifs, it->satellite, groupNameSatellite);
			}

			if(measureGroupsNameMap.find(groupName) != measureGroupsNameMap.end())
			{
				int id = measureGroupsNameMap[groupName];
				measureGroupIds[id].measureParamsSets.push_back(measureParamsSetsId++);
				it->measureGroupId = id;
			}
			else
			{
				measureGroupsNameMap[groupName] = measureGroupId++;
				it->measureGroupId = measureGroupsNameMap[groupName];
				MeasureGroup measureGroup;
				//TODO: next line adde by Sergeev need to check
				measureGroup.name = groupName;
				measureGroup.measureParamsSets.push_back(measureParamsSetsId++);
				measureGroupIds.push_back(measureGroup);
			}
			if(it->satellite != nullptr)
			{
				if(measureGroupsNameMap.find(groupNameSatellite) != measureGroupsNameMap.end())
				{
					int id = measureGroupsNameMap[groupNameSatellite];
					//measureGroupIds[id].measureParamsSets.push_back(measureParamsSetsId++);
					if (groupName != groupNameSatellite)
					{
						measureGroupIds[id].measureParamsSets.push_back(measureParamsSetsId);
					}
					it->satellite->measureGroupId = id;
				}
				else
				{
					measureGroupsNameMap[groupNameSatellite] = measureGroupId++;
					it->satellite->measureGroupId = measureGroupsNameMap[groupNameSatellite];
					MeasureGroup measureGroup;
					measureGroup.name = groupName;
					//measureGroup.measureParamsSets.push_back(measureParamsSetsId++);
					measureGroup.measureParamsSets.push_back(measureParamsSetsId);
					measureGroupIds.push_back(measureGroup);
				}
				//it->satellite->measureGroupId = it->measureGroupId;
			}
			++it;
		}
	}
	fs::ReadLineValue(ifs, currentId); // даже если 0, там будет 1 !
	currentId--;
}

class eq_binder: public binary_function<pair<string, int>, int, bool>
{
public:
	bool operator () (const pair<string, int> &param, int value) const
	{
		return param.second == value;
	}
};

void FMeasureParamsSet::Write(ofstream& ofs) const
{
	ofs << measureParamsSets.size() << std::endl;
	vector<MeasureParams>::const_iterator it = measureParamsSets.begin();
	while(it != measureParamsSets.end())
	{
		std::map<std::string, int>::const_iterator mit = 
			std::find_if( measureGroupsNameMap.begin(), measureGroupsNameMap.end(), 
			std::bind2nd(eq_binder(), it->measureGroupId));
		
		if( mit == measureGroupsNameMap.end())
		{
			exceptions::ThrowMessage("Measure param is outside of available Groups");
		}
		WriteMeasure(ofs, &(*it), mit->first);
		if(it->satellite != nullptr)
		{
			std::map<std::string, int>::const_iterator mit1 = 
				std::find_if( measureGroupsNameMap.begin(), measureGroupsNameMap.end(), 
				std::bind2nd(eq_binder(), it->satellite->measureGroupId));
		
			WriteMeasure(ofs, it->satellite, mit1->first);
		}
		++it;
	}
	ofs << this->currentId+1 << std::endl; // TODO: проверить для 0 и вообще этот факт
}

void FMeasureParamsSet::OutputGroup(ofstream& ofs, int id) const
{
	if (id < measureGroupIds.size())
	{
		vector<size_t>::const_iterator it = measureGroupIds[id].measureParamsSets.begin();
		map<int, vector<size_t> > outputClassToIds;

		while(it != measureGroupIds[id].measureParamsSets.end())
		{
			const MeasureParams& mp = measureParamsSets[*it];
			outputClassToIds[mp.outputClass].push_back(*it);
			++it;
		}

		map<int, vector<size_t> >::const_iterator mit = outputClassToIds.begin();
		while(mit != outputClassToIds.end())
		{
			it = mit->second.begin();
			while(it != mit->second.end())
			{
				const MeasureParams& mp = measureParamsSets[*it];
				ofs << mp.name << endl;
				ofs << setw(4) << mp.outputClass << ' ' << setw(2) << 
					mp.outputCode << ' ';

				if((mp.outputClass == 2) || (mp.outputClass == 3))
				{
					//			ofs << setw(4) << mp.subobject << ' ';
					ofs << setw(6) << mp.object << ' '; // номер соединительного элемента ?
					ofs << setw(4) << 1 << ' '; // любое число
				}
				else
				{
					ofs << setw(4) << mp.subobject << ' ';
					ofs << setw(6) << mp.object << ' '; // ??? nodeNumber and unknownVal
				}
				for(int i = 0; i < 6; i++)
					ofs << mp.dof[i] << ' ';
				for(int i = 0; i < 3; i++)
					ofs << mp.characteristic[i] << ' ';
				ofs << mp.isLogarifmic << ' '; // 14
				// TODO: разобраться с непонятными параметрами
				ofs << "0. 0. 0. "; // ???
				for(int i = 0; i < 3; i++)
					ofs << mp.params[i] << ' ';
				ofs << endl;
				++it;
			}
			++mit;
		}
	}
}

void FMeasureParamsSet::OutputGroup(ofstream& ofs, const string& groupName) const
{
	const auto& measureGroupIterator = measureGroupsNameMap.find(groupName);

	if (measureGroupIterator != measureGroupsNameMap.end())
	{
		int id = measureGroupIterator->second;

		OutputGroup(ofs, id);
	}
}

void FMeasureParamsSet::Output(ofstream& ofs) const
{
	OutputGroup(ofs, currentId);
}

void FMeasureParamsSet::BuildNameMap()
{
	measureGroupsNameMap.clear();
	int i=0;

	for (vector<MeasureGroup>::iterator it = measureGroupIds.begin(); it!=measureGroupIds.end(); ++it)
	{
		measureGroupsNameMap[it->name] = i;
		i++;
	}
}