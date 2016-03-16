#ifndef RESULTS_WRITER_H

#define RESULTS_WRITER_H


#include <string>
#include <fstream>
#include <memory>


using std::string;
using std::ofstream;


/** Класс для сброса результатов в файл
*
* @author Getmanskiy Victor, Sergeev Efim
*/
class ResultsWriter
{
public:

	ResultsWriter
		(
			const string& filename
		);

	~ResultsWriter();


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
	void WriteBuffer
		(
			const void* buffer,
			size_t bufferSize
		);

	/** Записать буфер double с преобразованием во float
	* @param buffer - буфер
	* @param bufferSize - размер буфера в элементах
	*/
	void WriteBufferToFloat
		(
			const double* buffer,
			size_t bufferSize
		);

private:
	ofstream _ofs;
	std::unique_ptr<float[]> _buf;
	void WriteFloatBuffer(const double* dbuf, float* fbuf, size_t bufferSize);

};


#endif // RESULTS_WRITER_H