#ifndef RLCControlReader_H
#define RLCControlReader_H
//#pragma GCC visibility push(hidden)
#include "ControlReader.h"

class RLCControlReader:public ControlReader
{
protected:
	
public:
	virtual OccRectilinearGrid* fillGrid(ifstream *fs);
	RLCControlReader(void);
	~RLCControlReader(void);
	/**
	* Функция читает бинарный файл со сжатым форматом хранения сетки
	* @param const char *fn имя файла
	* @return возвращет 1 если чтение произведено успешно
	*/
	//int ReadFile(const char *fn);
};
//#pragma GCC visibility pop
#endif
