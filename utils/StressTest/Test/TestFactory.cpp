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

SpecialSolvers::StressStrainStuff::SolverHandler SpecialSolversTest::StressStrainStuff::TestFactory::BuildBeam()
{
	stringstream str;
	str << "beam_stress_type_" << solverType << '_' << _gridParams._nx << 'x' << _gridParams._ny << 'x' << _gridParams._nz << '_' << ECodeToString(code);
	Uid(str.str());
	return MakeSolverBeam(
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

SpecialSolvers::StressStrainStuff::SolverHandler SpecialSolversTest::StressStrainStuff::TestFactory::BuildBeam2()
{
	stringstream str;
	str << "beam2_stress_type_" << solverType << '_' << _gridParams._nx << 'x' << _gridParams._ny << 'x' << _gridParams._nz << '_' << ECodeToString(code);
	Uid(str.str());
	return MakeSolverBeam2(
		_gridParams,
		_specialParams,
		_integrationParams,
		solverUid,
		faceSealed,
		faceSealed2,
		faceForced,
		force,
		forceDof,
		solverType
		);
}


SpecialSolvers::StressStrainStuff::SolverHandler SpecialSolversTest::StressStrainStuff::TestFactory::BuildPlate()
{
	stringstream str;
	str << "plate_stress_type_" << solverType << '_' << _gridParams._nx << 'x' << _gridParams._ny << 'x' << _gridParams._nz << '_' << ECodeToString(code);
	Uid(str.str());
	return MakePlateSolver(
		_gridParams,
		_specialParams,
		_integrationParams,
		solverUid,
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
		// TODO: убрать костыль
		//if (_gridParams._nx == 1) forceDof = dof_x;
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
	case xlrx:
		faceForced = SpecialSolversTest::face_left;
		faceSealed = SpecialSolversTest::face_right;
		std::swap(faceForced, faceSealed);
		_gridParams._nx = length;
		_gridParams._ny = sectionWidth;
		_gridParams._nz = sectionHeight;
		forceDof = dof_x;
		// TODO: убрать костыль
		//if (_gridParams._nx == 1) forceDof = dof_x;
		break;
	default:
		break;
	}
	return *this;
}

SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::Dims2(size_t length, size_t sectionWidth, size_t sectionHeight, ECode code)
{
	this->code = code;
	switch (code)
	{
	case xlr:
	case xrl:
		faceSealed = SpecialSolversTest::face_left;
		faceSealed2 = SpecialSolversTest::face_right;
		faceForced = SpecialSolversTest::face_front;
		_gridParams._nx = length;
		_gridParams._ny = sectionWidth;
		_gridParams._nz = sectionHeight;
		forceDof = dof_y;
		break;
	case yfb:
	case ybf:
		faceForced = SpecialSolversTest::face_left;
		faceSealed = SpecialSolversTest::face_front;
		faceSealed2 = SpecialSolversTest::face_back;
		_gridParams._nx = sectionWidth;
		_gridParams._ny = length;
		_gridParams._nz = sectionHeight;
		forceDof = dof_z;
		break;
	case ztb:
	case zbt:
		faceForced = SpecialSolversTest::face_front;
		faceSealed = SpecialSolversTest::face_top;
		faceSealed2 = SpecialSolversTest::face_bottom;
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


SpecialSolversTest::StressStrainStuff::TestFactory& SpecialSolversTest::StressStrainStuff::TestFactory::ForceDof(EDOF dof)
{
	forceDof = dof;
	return *this;
}
