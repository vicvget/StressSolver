#include "StressSolver.h"
#include "BlenderExporter.h"

#include "../Solvers/Stress/StressStrainSolverExports.h"
#include "../Solvers/Stress/StressStrainCppIterativeSolver.h"
#include "../FormatProviders/ResProvider/ProviderMpr.h"


#include <iostream>
#include "MprExporter.h"
#include "ChartsExporter.h"
#include "DummyExporter.h"
#include "FrameChartsExporter.h"


//#define NO_BLENDER
//#define NO_CHARTS
//#define NO_WRITE_RESULTS

#define DISABLE_OUTPUT

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
			//std::vector<size_t> ids = { 483, 484 };
			//std::shared_ptr<BaseExporter> chartsExporter = std::make_shared<ChartsExporter>(ssSolver, ids);

			std::shared_ptr<BaseExporter> blenderExporter = std::make_shared<DummyExporter>(ssSolver);
			std::shared_ptr<BaseExporter> chartsExporter = std::make_shared<DummyExporter>(ssSolver);
			std::shared_ptr<BaseExporter> frameChartsExporter = std::make_shared<DummyExporter>(ssSolver);
			//for (size_t id = 11 * 5; id < 11 * 6; id++)
			//std::vector<size_t> ids2;
			//for (size_t id = 1; id <= 10; id++)
			//{
			//	ids2.push_back(id);
			//}
			//std::shared_ptr<BaseExporter> frameChartsExporter = std::make_shared<FrameChartsExporter>(ssSolver, ids2);

			int iteration = 0;

			mprExporter->Init();
			blenderExporter->Init();
			chartsExporter->Init();
			frameChartsExporter->Init();

			int progress = 0;

			Stress::InitialSolve(hSolver);
			std::cout << "Initial solve completed" << std::endl << std::flush;

			for (int i = 0; i < integrationParams._nIterations; i++)
			{
				Stress::Solve(hSolver, integrationParams._nSubIterations);

				if (i*100. / integrationParams._nIterations > progress)
				{
					std::cout << std::endl << progress << "% " << std::flush;
					progress++;
					Stress::UpdateBufferWithOutput(hSolver);
				}
				else
				{
					Stress::UpdateBuffer(hSolver);
				}
				float time = (1 + iteration) * integrationParams._timeStep * integrationParams._nSubIterations;
				#ifndef DISABLE_OUTPUT
				mprExporter->WriteFrame(time);
				blenderExporter->WriteFrame(time);
				chartsExporter->WriteFrame(time);
				frameChartsExporter->WriteFrame(time);
				#endif
				iteration++;
			}
			mprExporter->WriteFrame(0);
			Stress::PrintTime(hSolver);

			mprExporter->Finalize();
			blenderExporter->Finalize();
			chartsExporter->Finalize();
			frameChartsExporter->Finalize();

		}	
	}
}
