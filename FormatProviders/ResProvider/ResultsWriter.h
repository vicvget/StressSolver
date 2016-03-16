#ifndef RESULTS_WRITER_H

#define RESULTS_WRITER_H


#include <string>
#include <fstream>
#include <memory>


using std::string;
using std::ofstream;


/** ����� ��� ������ ����������� � ����
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


	/** �������� ������ ��� ����������� ������
	* @param size - ������ ������
	*/
	void AllocateFloatBuffer
		(
			size_t size
		);

	/** ������� ����
	* @param filename - ��� �����
	*/
	void Open
		(
			const string& filename
		);

	/** ������� ����
	*/
	void Close();

	/** �������� �����
	* @param buffer - �����
	* @param bufferSize - ������ ������
	*/
	void WriteBuffer
		(
			const void* buffer,
			size_t bufferSize
		);

	/** �������� ����� double � ��������������� �� float
	* @param buffer - �����
	* @param bufferSize - ������ ������ � ���������
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