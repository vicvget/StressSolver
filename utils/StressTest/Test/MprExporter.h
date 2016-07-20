#pragma once 

#include "BaseExporter.h"
#include "../FormatProviders/ResProvider/ProviderMpr.h"
#include "CommonSolvers.h"

class MprExporter: public BaseExporter
{
	SpecialSolvers::IntegrationParams _integrationParams;
	ProviderMpr _writer;

public:
	MprExporter(Stress::StressStrainCppIterativeSolver* solver, SpecialSolvers::IntegrationParams integrationParams) :
		BaseExporter(solver), 
		_integrationParams(integrationParams) {};

#pragma region overriden
public:
	virtual void Init();
	virtual void WriteFrame(float time);
	virtual void Finalize();
#pragma endregion
};
