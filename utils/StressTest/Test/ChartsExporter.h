#include "BaseExporter.h"
#include <iostream>
#include <fstream>

class ChartsExporter: public BaseExporter
{
	ChartsExporter(Stress::StressStrainCppIterativeSolver* solver) : BaseExporter(solver) {};
	std::ofstream _ofs;
	virtual void Init()
	{
		_ofs.open("charts.txt");
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
	};

	size_t elementId = 2;

	virtual void WriteFrame(float time)
	{
		if (_state)
		_ofs << std::setprecision(5) << time << ' ' << std::setprecision(15)
			<< _solver->GetElementShift(elementId)[0] << ' '
			<< _solver->GetElementShift(elementId)[1] << ' '
			<< _solver->GetElementShift(elementId)[2] << ' '
			<< _solver->GetElementShiftAngular(elementId)[0] << ' '
			<< _solver->GetElementShiftAngular(elementId)[1] << ' '
			<< _solver->GetElementShiftAngular(elementId)[2] << ' '
			//<< _solver->GetElementVelocity(elementId)[0] << ' '
			//<< _solver->GetElementVelocity(elementId)[1] << ' '
			//<< _solver->GetElementVelocity(elementId)[2] << ' '
			//<< _solver->GetElementVelocityAngular(elementId)[0] << ' '
			//<< _solver->GetElementVelocityAngular(elementId)[1] << ' '
			//<< _solver->GetElementVelocityAngular(elementId)[2] << ' '
			<< _solver->GetElementAcceleration(elementId)[0] << ' '
			<< _solver->GetElementAcceleration(elementId)[1] << ' '
			<< _solver->GetElementAcceleration(elementId)[2] << ' '
			<< _solver->GetElementAccelerationAngular(elementId)[0] << ' '
			<< _solver->GetElementAccelerationAngular(elementId)[1] << ' '
			<< _solver->GetElementAccelerationAngular(elementId)[2] << ' '
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
	
	virtual void Finalize()
	{
		if (_ofs.is_open())
			_ofs.close();

	}


};