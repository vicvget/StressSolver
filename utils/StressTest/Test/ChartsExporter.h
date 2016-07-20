#pragma once 

#include "BaseExporter.h"
#include <iostream>
#include <fstream>

class ChartsExporter: public BaseExporter
{
	std::ofstream _ofs;
	size_t _elementId = 2;

public:
	ChartsExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};

#pragma region overriden
public:
	virtual void Init();
	virtual void WriteFrame(float time);
	virtual void Finalize();
#pragma endregion

};