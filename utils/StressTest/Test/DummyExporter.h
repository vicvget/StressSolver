#pragma once 

#include "BaseExporter.h"
#include <iostream>
#include <fstream>

class DummyExporter: public BaseExporter
{
public:
	DummyExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};

#pragma region overriden
public:
	virtual void Init();
	virtual void WriteFrame(float time);
	virtual void Finalize();
#pragma endregion

};