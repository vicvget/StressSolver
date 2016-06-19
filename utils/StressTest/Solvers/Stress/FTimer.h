#pragma once

//#define NOTIMER
//#undef _MSC_VER
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

	/* Subtract the ‘struct timeval’ values X and Y,
	storing the result in RESULT.
	Return 1 if the difference is negative, otherwise 0. */

	int	timeval_sub(struct timeval *result, struct timeval *x, struct timeval *y)
	{
		/* Perform the carry for the later subtraction by updating y. */
		if (x->tv_usec < y->tv_usec) {
			int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
			y->tv_usec -= 1000000 * nsec;
			y->tv_sec += nsec;
		}
		if (x->tv_usec - y->tv_usec > 1000000) {
			int nsec = (x->tv_usec - y->tv_usec) / 1000000;
			y->tv_usec += 1000000 * nsec;
			y->tv_sec -= nsec;
		}

		/* Compute the time remaining to wait.
		tv_usec is certainly positive. */
		result->tv_sec = x->tv_sec - y->tv_sec;
		result->tv_usec = x->tv_usec - y->tv_usec;

		/* Return 1 if result is negative. */
		return x->tv_sec < y->tv_sec;
	}

	int	timeval_add(struct timeval *result, struct timeval *x, struct timeval *y)
	{
		int nsec = 0;
		if (x->tv_usec + y->tv_usec > 1000000)
			nsec = (x->tv_usec + y->tv_usec) / 1000000;

		/* Compute the time remaining to wait.
		tv_usec is certainly positive. */
		result->tv_sec = nsec + x->tv_sec + y->tv_sec;
		result->tv_usec = x->tv_usec + y->tv_usec;

		/* Return 1 if result is negative. */
		return 0;
	}
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
		struct timeval current;
		gettimeofday(&current, NULL);

		struct timeval diff;
		timeval_sub(&diff, &current, &_current);
		timeval_add(&_total,&_total,&diff);
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