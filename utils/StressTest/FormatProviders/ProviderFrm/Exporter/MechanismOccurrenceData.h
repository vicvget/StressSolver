#ifndef MECHANISM_OCCURRENCE_DATA_H

#define MECHANISM_OCCURRENCE_DATA_H


#include <string>
#include <map>


namespace MechanismOccurrence
{

	using std::map;
	using std::string;


	// наименование тела
	typedef string BodyName;

	// идентификатор тела
	typedef int BodyId;

	// ???
	typedef map<BodyName, BodyId> MechanismOccurrenceData;


	/**
	* Загрузить данные ??? из файла с именем fileName
	* @param fileName - имя файла, из которого загружаются данные
	* @param data - загруженные данные
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool Load
		(
			const string& fileName,
			MechanismOccurrenceData& data
		);

}


#endif // MECHANISM_OCCURRENCE_DATA_H