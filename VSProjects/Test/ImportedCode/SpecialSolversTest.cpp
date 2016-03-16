#include "StressSolverTest.h"

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

int main1()
{
	using namespace SpecialSolversTest;
	StressStrainStuff::Test3x3x10(0);
	detect_memory_leaks(true);

	//StressStrainStuff::Test();
	//StressStrainStuff::Test3x3x10(3);
	//StressStrainStuff::Test3x3x10(2);
	//StressStrainStuff::Test3x1x1(-1);
	//StressStrainStuff::Test3x3x10(0);
	StressStrainStuff::Test3x3x10(0);
	//StressStrainStuff::Test1x3x10();
	//StressStrainStuff::Test3x3x10(2);
	//StressStrainStuff::TestSolveSystemOfLinearEquationsForStiffness();
	//ThermalStuff::Test();
	return 0;
}