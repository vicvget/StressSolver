#pragma once


#include <windows.h>
#include <string>


class SolverPerformanceCounter final
{
public:

	SolverPerformanceCounter();


	void Reset();


	double Time() const;

	void Print
		(
			const std::string& tag
		)	const;


private:

	LARGE_INTEGER _frequency;

	LARGE_INTEGER _current{};

};