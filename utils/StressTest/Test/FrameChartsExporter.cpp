#include "FrameChartsExporter.h"
#include <sstream>

void FrameChartsExporter::Init()
{
	if (_solver != nullptr)
	{
		SetInitialized();
	}
	else
	{
		SetError();
	}
}

void FrameChartsExporter::WriteFrame(float time)
{
	std::stringstream ss;
	ss << _solver->_uid << std::setprecision(5) << time << "_framechart.txt";
	//_ofs.
	_ofs.open(ss.str());
	if (_ofs.is_open() && _solver != nullptr)
	{
		_ofs << "x y z def\n";
		for (size_t i = 0; i < _elementIds.size(); i++)
		{
			size_t id = _elementIds[i];
			const double* coords = _solver->GetElementGridCoordinates(id);
			_ofs << std::setprecision(5) << coords[0] << ' ' << coords[1] << ' ' << coords[2] << ' ' <<
				std::setprecision(15) << _solver->GetElementDisplacement(id) << std::endl;
		}
		_ofs.close();
	}
}

void FrameChartsExporter::Finalize()
{
}
