//#ifdef USE_KNC
#pragma once

#include "StressStrainCppIterativeSolver.h"
#include "StressStrainCppIterativeSolverKNC.h"

#include <vector>
#include <immintrin.h>

using std::string;
using std::vector;

namespace Stress
{

	class StressStrainCppIterativeSolverKNC2
		:
		public StressStrainCppIterativeSolverKNC
	{
		//__m512d timeStep;
		//__m512d timeStep2;
		//__m512d timeStep4;
		//__m512d constantD2;
		//__m512d constantD6;
		const int regSize = 8;

	public:

		//virtual void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		// создает объект с заданными параметрами
		StressStrainCppIterativeSolverKNC2
			(
			double* params,
			int* links,
			int nLinks,
			double *nodes,
			int nNodes,
			double gridStep,
			double timeStep,
			int numThreads,
			int stride
			);

		virtual
			~StressStrainCppIterativeSolverKNC2();

#pragma region overriden


		virtual
		void CalculateStrains(size_t side, double* shiftStrains, double* velocityStrains, size_t nodeId1, size_t nodeId2) const;
		
		virtual
		void  CalculateForces();
		    
		void CrossProductTwice(double* v1, double* v2, double* v3, double* res) const;

#pragma endregion

		double df[12]; // debug
	};
}
//#endif