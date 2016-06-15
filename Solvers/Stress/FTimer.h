#pragma once


// include OS specific timing library
#ifdef _MSC_VER
// Windows
#include <Windows.h>
#else
// Linux
#include <sys/time.h>
#endif

#include <vector>
#include <iomanip>

class PerformanceCounter
{
private:
#ifdef _MSC_VER
	LARGE_INTEGER _freq;
	LARGE_INTEGER _current;
	LARGE_INTEGER _total;
#else
	// lower resolution but never wraps around:
	struct timeval _current;
	struct timeval _total;
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
		gettimeofday(&_current, NULL);
#endif
	}

	void Stop()
	{
#ifdef _MSC_VER
		LARGE_INTEGER current_time;
		QueryPerformanceCounter(&current_time);
		_total.QuadPart += (current_time.QuadPart - _current.QuadPart);
#else
		struct timeval _current;
		gettimeofday(&_current, NULL);
		_total += (current - _current);
#endif
	}

	double Get(bool isStop = false)
	{
		if (isStop) Stop();
		double current_time;

#ifdef _MSC_VER
		current_time = (double)(_total.QuadPart) / _freq.QuadPart;
#else
		current_time = (double)(_total.tv_sec + _total.tv_usec / 1e6);
#endif

		return current_time;
	}

	double Print(const char* tag, bool isStop = false)
	{
		double current_time = Get(isStop);
		std::cout << std::setw(10) << tag << std::setprecision(20) << current_time << std::endl;
		return current_time;
	}
};


class FTimer
{
	std::vector<PerformanceCounter>	_counters;
public:
	FTimer(){};
	void Add()
	{
		_counters.push_back(PerformanceCounter());
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