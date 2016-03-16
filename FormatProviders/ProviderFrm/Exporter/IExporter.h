#ifndef IEXPORTER_H

#define IEXPORTER_H


#include <string>


using std::string;


/** Интерфейс экспорта
*
* @author Getmanskiy Victor
*/
class IExporter
{
public:

	virtual ~IExporter(){};

	/** Экспорт по умолчанию
	*/
	virtual void Export() const = 0;

	/** Экспорт в заданный файл
	* @param file - файл
	*/
	virtual void ExportTo(const string& file) const = 0;

	/** Возвращает расширение
	* @return расширение
	*/
	virtual string GetExt() const = 0;

};


#endif // IEXPORTER_H
