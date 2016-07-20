#pragma once
#include <StressStrainCppIterativeSolver.h>

class BaseExporter
{
protected:
	bool _state;
	Stress::StressStrainCppIterativeSolver* _solver;
public:
	
	BaseExporter(Stress::StressStrainCppIterativeSolver* solver) : _solver(solver), _state(true) {};
	
	void SetInitialized() { _state = true; }
	void SetError() { _state = false; } 
	bool State() { return _state; }


	virtual void Init() = 0;
	virtual void WriteFrame(float time) = 0;
	virtual void Finalize() = 0;	
};
