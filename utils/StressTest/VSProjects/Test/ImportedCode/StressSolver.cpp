#include "StressSolver.h"

#include "../Stress/StressStrainSolverExports.h"
#include "../Stress/StressStrainCppIterativeSolver.h"
#include "../../FormatProviders/ResProvider/ProviderMpr.h"


#include <iostream>
#include "../BlenderExporter.h"


namespace Stress{
	class StressStrainCppIterativeSolver;
}

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

#ifndef NOBLENDER
			BlenderExporter exporter;			
			Stress::StressStrainCppIterativeSolver* ssSolver = static_cast<Stress::StressStrainCppIterativeSolver*>(hSolver);
			exporter.Init("blender.py", gridParams.NodesCount(), ssSolver->vecStride, ssSolver->GetElementShift(0), gridParams._gridStep);
			exporter.WriteHeader();
			exporter.WriteBody();
#endif

			Stress::UpdateBuffer(hSolver, scaleResults);
			float* data = Stress::GetMemoryPointer(hSolver);
			int dataSize = Stress::GetMemorySize(hSolver);
			writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep * integrationParams._nSubIterations);
			int cP = 0;

			ofstream dofs("charts.txt");
			dofs << "t x y z rx ry rz vx vy vz wx wy wz ax ay az ex ey ez\n";
			for (int i = 0; i < integrationParams._nIterations; i++)
			{
				Stress::Solve(hSolver, integrationParams._nSubIterations);
#ifndef NO_WRITE_RESULTS
				//				if (i % 100 == 0)
				{
					Stress::UpdateBuffer(hSolver, scaleResults);
					//float y = Stress::GetData(hSolver, nx - 1, 0, 0);
					//float stress = Stress::GetStressData(hSolver, nx - 1);
					//file << x << " " << y << " " << stress << std::endl;
#ifndef NOBLENDER
					exporter.AddFrame(ssSolver->GetElementShift(0), ssSolver->GetDataRotaionMtx());
#endif
					writer.WriteFrame(data, dataSize, iteration * integrationParams._timeStep * integrationParams._nSubIterations);

					int elementId = 2;
					dofs << i << ' '
						<< ssSolver->GetElementShift(elementId)[0] << ' '
						<< ssSolver->GetElementShift(elementId)[1] << ' '
						<< ssSolver->GetElementShift(elementId)[2] << ' '
						<< ssSolver->GetElementShiftAngular(elementId)[0] << ' '
						<< ssSolver->GetElementShiftAngular(elementId)[1] << ' '
						<< ssSolver->GetElementShiftAngular(elementId)[2] << ' '
						<< ssSolver->GetElementVelocity(elementId)[0] << ' '
						<< ssSolver->GetElementVelocity(elementId)[1] << ' '
						<< ssSolver->GetElementVelocity(elementId)[2] << ' '
						<< ssSolver->GetElementVelocityAngular(elementId)[0] << ' '
						<< ssSolver->GetElementVelocityAngular(elementId)[1] << ' '
						<< ssSolver->GetElementVelocityAngular(elementId)[2] << ' '
						<< ssSolver->GetElementAcceleration(elementId)[0] << ' '
						<< ssSolver->GetElementAcceleration(elementId)[1] << ' '
						<< ssSolver->GetElementAcceleration(elementId)[2] << ' '
						<< ssSolver->GetElementAccelerationAngular(elementId)[0] << ' '
						<< ssSolver->GetElementAccelerationAngular(elementId)[1] << ' '
						<< ssSolver->GetElementAccelerationAngular(elementId)[2] << '\n';
					
				}

#endif
				if (i*100. / integrationParams._nIterations > cP)
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
#ifndef NOBLENDER
			exporter.WriteFooter();
#endif
		}	
	}

}
