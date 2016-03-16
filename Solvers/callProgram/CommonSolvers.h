#pragma once


namespace SpecialSolvers
{

	struct IntegrationParams
	{

		int _nIterations;

		int _nSubIterations;

		float _initialTime;

		float _timeStep;


		IntegrationParams();

		float TotalTime() const;

	};


	struct GridParams
	{

		double _gridStep;

		int _nx;

		int _ny;

		int _nz;

		int _gridNodesCount{};


		GridParams();

		int NodesCount() const;

	};

}