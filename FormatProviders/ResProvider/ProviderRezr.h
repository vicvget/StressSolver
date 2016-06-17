#ifndef PROVIDER_REZR_H

#define PROVIDER_REZR_H


#include "FileStream.h"
#include "AbstractProviderRezr.h"


/**
* Провайдер файла результатов
*/
class ProviderRezr
	:
		public AbstractProviderRezr
{
public:

	// Конструкторы/деструктор

	/**
	* Конструктор
	* @param rezrName - имя файла с результатами
	* @param fadresName - имя файла со структурой результатов
	*/
	ProviderRezr
		(
			const string& rezrName,
			const string& fadresName
		);
	
	/**
	* Деструктор
	*/
	virtual
		~ProviderRezr();// override;


	// Работа с файлом результатов

	/**
	* Открытие файла результатов
	* @return false в случае ошибки
	*/
	virtual
	bool OpenStream() override;

	/**
	* Закрыть поток результатов
	*/
	virtual
	void CloseStream() override;

	/**
	* Считывание заголовка
	*/
	virtual
	void ReadHeader() override;

	/**
	* Чтение кадра
	* @return false, если кадр не прочитан целиком или конец потока
	*/
	virtual
	bool TryReadFrame() override;


protected:

	// файловый поток для чтения файла результатов
	PowerStationFileStream* _powerStationFileStream;

};


#endif // PROVIDER_REZR_H