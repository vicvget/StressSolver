#ifndef RESULTS_RAW_READER_H

#define RESULTS_RAW_READER_H


#include <string>
#include <fstream>
#include <memory>


using std::string;
using std::ifstream;
using std::unique_ptr;

/** Класс для чтения результатов
*
* @author Getmanskiy Victor
*/
class ResultsRawReader
{
public:

	ResultsRawReader
		(
			const string& filename
		);

	~ResultsRawReader();


	/** Выделить память для внутреннего буфера
	* @param size - размер буфера
	*/
	void AllocateFloatBuffer
		(
			size_t size
		);

	/** Открыть файл
	* @param filename - имя файла
	*/
	void Open
		(
			const string& filename
		);

	/** Закрыть файл
	*/
	void Close();

	/** Записать буфер
	* @param buffer - буфер
	* @param bufferSize - размер буфера
	*/
	bool ReadBuffer
		(
			void* buffer,
			size_t bufferSize
		);

private:
	ifstream _ifs;
	unique_ptr<float> _buf;
};


#endif // RESULTS_RawReader_H