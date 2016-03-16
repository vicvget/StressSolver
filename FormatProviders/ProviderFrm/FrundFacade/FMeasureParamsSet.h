#ifndef FMEASURE_PARAMS_SET_H

#define FMEASURE_PARAMS_SET_H


#include "FParamsSet.h"

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>


using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::map;


struct MeasureParams
{
	
	MeasureParams();

	MeasureParams
		(
			const MeasureParams* mp
		);

	MeasureParams
		(
			const MeasureParams& mp
		);

	~MeasureParams();

	MeasureParams& operator = (const MeasureParams& mp);

	void Copy(const MeasureParams& mp);

	void CopyMembers(const MeasureParams& mp);


	MeasureParams* satellite;
	string name;
	size_t object;		//body, link
	size_t subobject;	//for body - node, other = 12
	size_t outputClass;	// 1-8 1-move, 2-forces, 3-deformations
	size_t outputCode;	// 1-30?
	char dof[6];
	size_t characteristic[3]; // for body 0,1 - displacement, velocity, acceleration, for force spring, damper and result
	bool isLogarifmic;		  // logarifmic axis
	float params[3];		  // params of chart
	size_t object1;
	size_t subobject1;
	size_t windowType; // 1-5 
	size_t fftBlockLength; // 2^n
	bool isMovingCS;
	int measureGroupId;

};

struct MeasureGroup
{
	string name;
	vector<size_t> measureParamsSets;
};

class FMeasureParamsSet : public FParamsSet
{
	vector<MeasureParams> measureParamsSets;
	vector<MeasureGroup> measureGroupIds;
	map<string, int> measureGroupsNameMap;
public:
	FMeasureParamsSet();
	~FMeasureParamsSet();

	size_t CurrentId() const { return currentId; }
	void CurrentId(size_t val) { currentId = val; }

	size_t Size() const {return measureGroupIds.size();}
	void Read(ifstream& ofs);
	void Write(ofstream& ofs) const;
	void BuildNameMap();
	/**
	* Генерирует файлы, необходимые для постобработки: uprf, ...
	*/
	void Output(ofstream& ofs) const;
	void OutputGroup(ofstream& ofs, const string& groupName = "common") const;
	void OutputGroup(ofstream& ofs, int id) const;


	vector<MeasureParams> MeasureParamsSets() const { return measureParamsSets; }
	void MeasureParamsSets(vector<MeasureParams> val) { measureParamsSets = val; }

	vector<MeasureGroup> MeasureGroupIds() const { return measureGroupIds; }
	void MeasureGroupIds(vector<MeasureGroup> val) { measureGroupIds = val; }
};


#endif // FMEASURE_PARAMS_SET_H