#pragma once


#include <windows.h>
#include <vector>


class PerformanceCounter
{
private:
	LARGE_INTEGER _freq;
	LARGE_INTEGER _current;
	LARGE_INTEGER _total;
public: 
	PerformanceCounter()
	{
		QueryPerformanceFrequency(&_freq);
		_total.QuadPart = 0LL;
	}
	void Start()
	{
		QueryPerformanceCounter(&_current);
	}

	void Stop()
	{
		LARGE_INTEGER current_time;	
		QueryPerformanceCounter(&current_time);
		_total.QuadPart+=(current_time.QuadPart - _current.QuadPart);
	}

	double Get(bool isStop = false)
	{
		if(isStop) Stop();
		double current_time = (double)(_total.QuadPart)/_freq.QuadPart;
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
		if(id < _counters.size())
			_counters[id].Start();
	}

	void Stop(unsigned int id)
	{
		if(id < _counters.size())
			_counters[id].Stop();
	}

	double Print(unsigned int id, const char* tag)
	{
		return id < _counters.size() ? _counters[id].Print(tag) : 0.;
	}
};