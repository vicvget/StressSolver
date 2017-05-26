#define SQR(x) ((x) * (x))
#define MAX(x, y) ((x) > (y) ? x : y)

#define MeasuredRun(TIMER, COMMAND) \
	_testTimer.Start(TIMER); \
	COMMAND; \
	_testTimer.Stop(TIMER);

#define ALIGNMENT 64 // KNC