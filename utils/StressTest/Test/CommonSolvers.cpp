#include "CommonSolvers.h"


namespace SpecialSolvers
{

	// struct IntegrationParams

	IntegrationParams::IntegrationParams()
		:
			_nIterations(100),
			_nSubIterations(1),
			_initialTime(0.0f),
			_timeStep(1e-2f)
	{
	}

	float IntegrationParams::TotalTime() const
	{
		return _nIterations * _nSubIterations * _timeStep;
	}


	// struct GridParams

	GridParams::GridParams()
		:
			_gridStep(0.01),
			_nx(30),
			_ny(3),
			_nz(1)
	{
	}

	int GridParams::NodesCount() const
	{
		return _nx * _ny * _nz;
	}

}