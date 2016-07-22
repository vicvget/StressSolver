#include "ChartsExporter.h"

void ChartsExporter::Init()
{
	_ofs.open(_solver->_uid + "_charts.txt");
	if (_ofs.is_open() && _solver != nullptr)
	{
		SetInitialized();
		//	dofs << "t x y z rx ry rz vx vy vz wx wy wz ax ay az ex ey ez dx dy dz drx dry drz fx fy fz frx fry frz\n";
		_ofs << "t x y z rx ry rz ax ay az ex ey ez\n";
	}
	else
	{
		SetError();
	}
}

void ChartsExporter::WriteFrame(float time)
{
	if (_state)
		_ofs << std::setprecision(5) << time << ' ' << std::setprecision(15)
			<< _solver->GetElementShift(_elementId)[0] << ' '
			<< _solver->GetElementShift(_elementId)[1] << ' '
			<< _solver->GetElementShift(_elementId)[2] << ' '
			<< _solver->GetElementShiftAngular(_elementId)[0] << ' '
			<< _solver->GetElementShiftAngular(_elementId)[1] << ' '
			<< _solver->GetElementShiftAngular(_elementId)[2] << ' '
			//<< _solver->GetElementVelocity(_elementId)[0] << ' '
			//<< _solver->GetElementVelocity(_elementId)[1] << ' '
			//<< _solver->GetElementVelocity(_elementId)[2] << ' '
			//<< _solver->GetElementVelocityAngular(_elementId)[0] << ' '
			//<< _solver->GetElementVelocityAngular(_elementId)[1] << ' '
			//<< _solver->GetElementVelocityAngular(_elementId)[2] << ' '
			<< _solver->GetElementAcceleration(_elementId)[0] << ' '
			<< _solver->GetElementAcceleration(_elementId)[1] << ' '
			<< _solver->GetElementAcceleration(_elementId)[2] << ' '
			<< _solver->GetElementAccelerationAngular(_elementId)[0] << ' '
			<< _solver->GetElementAccelerationAngular(_elementId)[1] << ' '
			<< _solver->GetElementAccelerationAngular(_elementId)[2] << ' '
			//<< ssSolver->df[0] << ' '
			//<< ssSolver->df[1] << ' '
			//<< ssSolver->df[2] << ' '
			//<< ssSolver->df[3] << ' '
			//<< ssSolver->df[4] << ' '
			//<< ssSolver->df[5] << ' '
			//<< ssSolver->df[6] << ' '
			//<< ssSolver->df[7] << ' '
			//<< ssSolver->df[8] << ' '
			//<< ssSolver->df[9] << ' '
			//<< ssSolver->df[10] << ' '
			//<< ssSolver->df[11]
			<< std::endl;
}

void ChartsExporter::Finalize()
{
	if (_ofs.is_open())
		_ofs.close();

}
