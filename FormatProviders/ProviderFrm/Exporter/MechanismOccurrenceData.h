#ifndef MECHANISM_OCCURRENCE_DATA_H

#define MECHANISM_OCCURRENCE_DATA_H


#include <string>
#include <map>


namespace MechanismOccurrence
{

	using std::map;
	using std::string;


	// ������������ ����
	typedef string BodyName;

	// ������������� ����
	typedef int BodyId;

	// ???
	typedef map<BodyName, BodyId> MechanismOccurrenceData;


	/**
	* ��������� ������ ??? �� ����� � ������ fileName
	* @param fileName - ��� �����, �� �������� ����������� ������
	* @param data - ����������� ������
	* @return ������� �������� (true) ��� ���������� (false) ��������
	*/
	bool Load
		(
			const string& fileName,
			MechanismOccurrenceData& data
		);

}


#endif // MECHANISM_OCCURRENCE_DATA_H