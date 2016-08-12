#include "StressSolver.h"
#include "BlenderExporter.h"

#include "../Solvers/Stress/StressStrainSolverExports.h"
#include "../Solvers/Stress/StressStrainCppIterativeSolver.h"
#include "../FormatProviders/ResProvider/ProviderMpr.h"


#include <iostream>
#include "MprExporter.h"
#include "ChartsExporter.h"
#include "DummyExporter.h"


//#define NO_BLENDER
//#define NO_CHARTS
//#define NO_WRITE_RESULTS

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
				SolverHandler hSolver,
//				const string& fileResults,
//				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			)
		{
			Stress::StressStrainCppIterativeSolver* ssSolver = static_cast<Stress::StressStrainCppIterativeSolver*>(hSolver);

			std::shared_ptr<BaseExporter> mprExporter = std::make_shared<MprExporter>(ssSolver, integrationParams);
			//std::shared_ptr<BaseExporter> blenderExporter = std::make_shared<BlenderExporter>(ssSolver);
			//std::shared_ptr<BaseExporter> chartsExporter = std::make_shared<ChartsExporter>(ssSolver);

			std::shared_ptr<BaseExporter> blenderExporter = std::make_shared<DummyExporter>(ssSolver);
			std::shared_ptr<BaseExporter> chartsExporter = std::make_shared<DummyExporter>(ssSolver);


			int iteration = 0;

			mprExporter->Init();
			blenderExporter->Init();
			chartsExporter->Init();

			int progress = 0;

			Stress::InitialSolve(hSolver);
			std::cout << "Initial solve completed" << std::endl << std::flush;

			for (int i = 0; i < integrationParams._nIterations; i++)
			{
				Stress::Solve(hSolver, integrationParams._nSubIterations);
				Stress::UpdateBuffer(hSolver);
				double time = (1 + iteration) * integrationParams._timeStep * integrationParams._nSubIterations;

				mprExporter->WriteFrame(time);
				blenderExporter->WriteFrame(time);
				chartsExporter->WriteFrame(time);

				if (i*100. / integrationParams._nIterations > progress)
				{
					std::cout << progress << '%' << std::endl << std::flush;
					progress++;
				}
				iteration++;
			}

			mprExporter->Finalize();
			blenderExporter->Finalize();
			chartsExporter->Finalize();

		}	
	}

}
