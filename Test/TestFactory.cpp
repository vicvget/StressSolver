#include "TestFactory.h"

string SpecialSolversTest::StressStrainStuff::ECodeToString(ECode code)
{
	switch (code)
	{
	case xlr: return "xlr";
	case xrl: return "xrl";
	case yfb: return "yfb";
	case ybf: return "ybf";
	case ztb: return "ztb";
	case zbt: return "zbt";
	}
	return "";
}

SpecialSolvers::GridParams SpecialSolversTest::StressStrainStuff::TestFactory::GridParams()
{ return _gridParams; }

SpecialSolvers::IntegrationParams SpecialSolversTest::StressStrainStuff::TestFactory::IntegrationParams()
{ return _integrationParams; }

SpecialSolvers::StressStrainStuff::SpecialParams SpecialSolversTest::StressStrainStuff::TestFactory::SpecialParams()
{ return _specialParams; }

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::E(float val)
{
	_specialParams._E = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Density(float val)
{
	_specialParams._density = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Damping(float val)
{
	_specialParams._dampingRatio = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::ScaleFactor(float val)
{
	_specialParams._scaleFactor = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::IterationsCount(size_t val)
{
	_integrationParams._nIterations = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::SubIterationsCount(size_t val)
{
	_integrationParams._nSubIterations = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::TimeStep(float val)
{
	_integrationParams._timeStep = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::GridStep(float val)
{
	_gridParams._gridStep = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Force(float val)
{
	force = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Uid(const string& val)
{
	solverUid = val;
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::SolverType(int val)
{
	solverType = val;
	return *this;
}

SpecialSolvers::StressStrainStuff::SolverHandler SpecialSolversTest::StressStrainStuff::TestFactory::Build()
{
	stringstream str;
	str << "test_stress_type_" << solverType << '_' << ECodeToString(code);
	Uid(str.str());
	return MakeSolver(
		_gridParams,
		_specialParams,
		_integrationParams,
		solverUid,
		faceSealed,
		faceForced,
		force,
		forceDof,
		solverType
	);
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Dims(size_t length, size_t sectionWidth, size_t sectionHeight, ECode code)
{
	this->code = code;
	switch (code)
	{
	case xlr:
		faceForced = SpecialSolversTest::face_left;
		faceSealed = SpecialSolversTest::face_right;
	case xrl:
		std::swap(faceForced, faceSealed);
		_gridParams._nx = length;
		_gridParams._ny = sectionWidth;
		_gridParams._nz = sectionHeight;
		forceDof = dof_y;
		break;
	case yfb:
		faceForced = SpecialSolversTest::face_front;
		faceSealed = SpecialSolversTest::face_back;
	case ybf:
		std::swap(faceForced, faceSealed);
		_gridParams._nx = sectionWidth;
		_gridParams._ny = length;
		_gridParams._nz = sectionHeight;
		forceDof = dof_z;
		break;
	case ztb:
		faceForced = SpecialSolversTest::face_top;
		faceSealed = SpecialSolversTest::face_bottom;
	case zbt:
		std::swap(faceForced, faceSealed);
		_gridParams._nx = sectionHeight;
		_gridParams._ny = sectionWidth;
		_gridParams._nz = length;
		forceDof = dof_x;
		break;
	default:
		break;
	}
	return *this;
}
