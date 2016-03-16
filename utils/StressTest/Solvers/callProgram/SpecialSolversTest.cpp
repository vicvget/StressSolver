#include "OilMistSolverTest.h"
#include "StressSolverTest.h"
#include "ThermalSolverTest.h"

#include <crtdbg.h>


void detect_memory_leaks(bool on_off)
{
	int flags = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);

	if (!on_off)
		flags &= ~_CRTDBG_LEAK_CHECK_DF;
	else
	{
		flags |= _CRTDBG_LEAK_CHECK_DF;
		_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
		_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	}

	_CrtSetDbgFlag(flags);
}

int main()
{
	using namespace SpecialSolversTest;

	detect_memory_leaks(true);
	//OilMistStuff::Test();

	StressStrainStuff::Test();

	//StressStrainStuff::TestSolveSystemOfLinearEquationsForStiffness();
	//ThermalStuff::Test();
}