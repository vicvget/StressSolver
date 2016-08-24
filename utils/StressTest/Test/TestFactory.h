#pragma once 

#include <sstream>
#include <vector>

#include "CommonSolversTest.h"
#include "StressSolverTest.h"

using std::string;
using std::vector;
using std::stringstream;

namespace SpecialSolversTest
{
	namespace StressStrainStuff
	{
			enum ECode
			{
				xlr,
				xrl,
				yfb,
				ybf,
				ztb,
				zbt,
				xlrx
			};

		string ECodeToString(ECode code);

		class TestFactory
		{
			SpecialParams _specialParams;
			IntegrationParams _integrationParams;
			GridParams _gridParams;
			ECode code;

		public:

			GridParams GridParams();
			IntegrationParams IntegrationParams();
			SpecialParams SpecialParams();

			SolverHandler BuildBeam();
			SolverHandler BuildBeam2();
			SolverHandler BuildPlate();

#pragma region setters
			TestFactory& E(float val);

			TestFactory& Density(float val);

			TestFactory& Damping(float val);

			TestFactory& ScaleFactor(float val);

			TestFactory& IterationsCount(size_t val);

			TestFactory& SubIterationsCount(size_t val);

			TestFactory& TimeStep(float val);

			TestFactory& GridStep(float val);

			TestFactory& Force(float val);

			TestFactory& Uid(const string& val);

			TestFactory& SolverType(int val);

#pragma endregion

			EDOF forceDof;
			EFACE faceSealed;
			EFACE faceSealed2;
			EFACE faceForced;
			float force;
			string solverUid;
			int solverType;

			TestFactory& Dims(size_t length, size_t sectionWidth, size_t sectionHeight, ECode code);
			TestFactory& Dims2(size_t length, size_t sectionWidth, size_t sectionHeight, ECode code);
			TestFactory& ForceDof(EDOF dof);
		};
	} // namespace StressStrainStuff
} // namespace SpecialSolversTest
