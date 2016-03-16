#include "PerformanceCounter.h"

#include <iostream>
#include <iomanip>


// class SolverPerformanceCounter

SolverPerformanceCounter::SolverPerformanceCounter()
{
	QueryPerformanceFrequency(&_frequency);
}

void SolverPerformanceCounter::Reset()
{
	QueryPerformanceCounter(&_current);
}

double SolverPerformanceCounter::Time() const
{
	LARGE_INTEGER time;

	QueryPerformanceCounter(&time);

	return double(time.QuadPart - _current.QuadPart) / _frequency.QuadPart;
}

void SolverPerformanceCounter::Print
	(
		const std::string& tag
	)	const
{
	std::cout << std::setw(10) << tag << std::setprecision(20) << Time() << std::endl;
}