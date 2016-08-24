#pragma once 

#include "BaseExporter.h"
#include <iostream>
#include <fstream>
#include <vector>

class FrameChartsExporter: public BaseExporter
{
	std::ofstream _ofs;
	vector<size_t> _elementIds;

public:
	FrameChartsExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};
	void AddId(size_t id) { _elementIds.push_back(id); }
#pragma region overriden
public:
	virtual void Init();
	virtual void WriteFrame(float time);
	virtual void Finalize();
#pragma endregion

};