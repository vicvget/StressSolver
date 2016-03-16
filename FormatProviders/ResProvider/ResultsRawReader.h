#ifndef RESULTS_RAW_READER_H

#define RESULTS_RAW_READER_H


#include <string>
#include <fstream>
#include <memory>


using std::string;
using std::ifstream;
using std::unique_ptr;

/** ����� ��� ������ �����������
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