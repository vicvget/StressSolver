#include "ThermalSolver.h"

#include "../Thermal_new/ThermalSolverExports.h"
#include "../../FormatProviders/ResProvider/ProviderMpr.h"

#include <iostream>


namespace SpecialSolvers
{

	namespace ThermalStuff
	{

		// struct SpecialParams

		SpecialParams::SpecialParams()
			:
				_tConduction(46.5),
				_tCapacity(460),
				_density(7800),
				_tBody(20)
		{
		}

		SpecialParams::SpecialParams
			(
				double tConduction,
				double tCapacity,
				double density,
				double tBody
			)
			:
				_tConduction(tConduction),
				_tCapacity(tCapacity),
				_density(density),
				_tBody(tBody)
		{
		}

		void SpecialParams::GetParams
			(
				double* params
			)	const
		{
			params[0] = _tConduction;
			params[1] = _tCapacity;
			params[2] = _density;
			params[3] = _tBody;
		}


		void Solve
			(
				ThermalSolver hSolver,
				const string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			)
		{
			MultiphysicsResultsHeader mprHeader
				(
					SolverTypes::ST_Thermal,
					gridParams.NodesCount(),
					integrationParams._initialTime,
					integrationParams.TotalTime()
				);
			ProviderMpr writer;
			writer.InitWriter(fileResults, &mprHeader);
			int iteration = 0;

			Thermal::UpdateReturnedBuffer(hSolver);
			const float* data = Thermal::GetReturnedBuffer(hSolver);
			int dataSize = Thermal::GetReturnedBufferSize(hSolver);
			writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep);
			int cP = 0;
			for (int i = 0; i < integrationParams._nIterations; i++)
			{
				Thermal::Solve(hSolver);
				if ((i + 1) % integrationParams._nSubIterations == 0)
				{
					Thermal::UpdateReturnedBuffer(hSolver);
					writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep);
				}
				if (i * 100. / integrationParams._nIterations > cP)
				{
					std::cout << cP << '%' << std::endl;
					cP++;
				}
				iteration++;
			}
		}
	}

}