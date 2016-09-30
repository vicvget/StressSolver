#pragma once 

#include "BaseExporter.h"
#include <iostream>
#include <fstream>

class ChartsExporter: public BaseExporter
{
	std::ofstream _ofs;
	vector<size_t> _elementIds;
	bool _firstTime;

public:
	ChartsExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};
	ChartsExporter(Stress::StressStrainCppIterativeSolver* solver, const vector<size_t>& ids) : BaseExporter(solver), _elementIds(ids) {};

#pragma region overriden
public:
	virtual void Init();
	virtual void WriteFrame(float time);
	virtual void Finalize();
#pragma endregion

};