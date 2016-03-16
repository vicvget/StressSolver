#include "StressSolver.h"

#include "../Stress/StressStrainSolverExports.h"
#include "../../FormatProviders/ResProvider/ProviderMpr.h"

#include <iostream>


namespace SpecialSolvers
{

	namespace StressStrainStuff
	{

		// struct SpecialParams

		SpecialParams::SpecialParams()
			:
				_E(2.07e11),
				_dampingRatio(1e-1),
				_density(7800),
				_scaleFactor(1e8)
		{
		}

		SpecialParams::SpecialParams
			(
				double E,
				double dampingRatio,
				double density,
				double scaleFactor
			)
			:
				_E(E),
				_dampingRatio(dampingRatio),
				_density(density),
				_scaleFactor(scaleFactor)
		{
		}

		void SpecialParams::GetParams
			(
				double* params
			)	const
		{
			params[0] = _E;
			params[1] = _dampingRatio;
			params[2] = _density;
			params[3] = _scaleFactor;
		}


		void Solve
			(
				StressStrainSolver hSolver,
				const string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			)
		{
			MultiphysicsResultsHeader mprHeader
				(
					SolverTypes::ST_StressStrain,
					gridParams.NodesCount(),
					integrationParams._initialTime,
					integrationParams.TotalTime()
				);

			ProviderMpr writer;
			writer.InitWriter(fileResults, &mprHeader);
			int iteration = 0;
			static const double scaleResults = 1.0;

			Stress::UpdateBuffer(hSolver, scaleResults);
			float* data = Stress::GetMemoryPointer(hSolver);
			int dataSize = Stress::GetMemorySize(hSolver);
			writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep * integrationParams._nSubIterations);
			int cP = 0;
			for (int i = 0; i < integrationParams._nIterations; i++)
			{
				Stress::Solve(hSolver, integrationParams._nSubIterations);
//				if (i % 100 == 0)
				{
					Stress::UpdateBuffer(hSolver, scaleResults);
					//float y = Stress::GetData(hSolver, nx - 1, 0, 0);
					//float stress = Stress::GetStressData(hSolver, nx - 1);
					//file << x << " " << y << " " << stress << std::endl;
					writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep * integrationParams._nSubIterations);
				}
				if(i*100./integrationParams._nIterations > cP)
				{
					std::cout << cP << '%' << std::endl;
					cP++;
				}
				//if(curForceIteration<forceIterations)
				//{
				//	curForceIteration++;
				//}
				//else
				//{
				//	Stress::ChangeBoundaryParam(hSolver, 1, 2, 0);
				//	//Stress::ChangeBoundaryParam(hSolver, 1, 0, 0);
				//}
				//if(i==(int)(nIterations/2))
				//{
				//	curForceIteration=0;
				//	Stress::ChangeBoundaryParam(hSolver,1,2,-1e3);
				//	Stress::ChangeBoundaryParam(hSolver,1,0,1e3);
				//}
				iteration++;
			}
		}	
	}

}