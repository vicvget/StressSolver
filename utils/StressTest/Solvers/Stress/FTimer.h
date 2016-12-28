#pragma once

//#define NOTIMER
//#undef _MSC_VER
//struct timeval
//{
//	long tv_sec;
//	long tv_usec;
//};
//
//void gettimeofday(timeval*, void*)
//{
//	
//}

// include OS specific timing library
#ifdef _MSC_VER
// Windows
#include <Windows.h>
#else
// Linux
#include <ctime>
#endif

#include <vector>
#include <iomanip>
#include <iostream>
#include <omp.h>

class PerformanceCounter
{
private:
#ifdef _MSC_VER
	LARGE_INTEGER _freq;
	LARGE_INTEGER _current;
	LARGE_INTEGER _total;
#else
	// lower resolution but never wraps around:
	clock_t _current;
	long _total;
#endif


public:
	PerformanceCounter()
	{
#ifdef _MSC_VER
		QueryPerformanceFrequency(&_freq);
		_total.QuadPart = 0LL;
#endif
	}
	void Start()
	{
#ifdef _MSC_VER
		QueryPerformanceCounter(&_current);
#else
		_total = 0;
#endif
	}

	void Stop()
	{
#ifdef _MSC_VER
		LARGE_INTEGER current_time;
		QueryPerformanceCounter(&current_time);
		_total.QuadPart += (current_time.QuadPart - _current.QuadPart);
#else
		_total += clock() - _current;
#endif
	}

	double Get(bool isStop = false)
	{
		if (isStop) Stop();
		double current_time;

#ifdef _MSC_VER
		current_time = (double)(_total.QuadPart) / _freq.QuadPart;
#else
		current_time = (double)_total / CLOCKS_PER_SEC;
#endif

		return current_time;
	}

	double Print(const char* tag, bool isStop = false)
	{
		double current_time = Get(isStop);
		std::cout << std::setw(15) << tag << std::setprecision(15) << current_time << std::endl;
		return current_time;
	}
};

class PerformanceOmpCounter
{
	double _current;
	double _total;
public:
	PerformanceOmpCounter()
	{
		_total = 0;
	}

	void Start()
	{
		_current = omp_get_wtime();
	}

	double Get(bool isStop = false)
	{
		if (isStop) Stop();
		double current_time = (double)_total;
		return current_time;
	}

	void Stop()
	{
		_total += omp_get_wtime() -_current;
	}

	double Print(const char* tag, bool isStop = false)
	{
		double current_time = Get(isStop);
		std::cout << std::setw(15) << tag << std::setprecision(15) << current_time << std::endl;
		return current_time;
	}
};



class FTimer
{
	std::vector<PerformanceOmpCounter>	_counters;
public:
	FTimer(){};
	void Add()
	{
		_counters.push_back(PerformanceOmpCounter());
	}
	void Allocate(unsigned int num)
	{
		_counters.resize(num);
	}

	double Get(unsigned int id)
	{
		return id > _counters.size() ? _counters[id].Get() : -1;
	}

	void Start(unsigned int id)
	{
		if (id < _counters.size())
			_counters[id].Start();
	}

	void Stop(unsigned int id)
	{
		if (id < _counters.size())
			_counters[id].Stop();
	}

	double Print(unsigned int id, const char* tag)
	{
		return id < _counters.size() ? _counters[id].Print(tag) : 0.;
	}
};