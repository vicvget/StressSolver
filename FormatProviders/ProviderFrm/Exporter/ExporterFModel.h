#ifndef EXPORTER_F_MODEL_H

#define EXPORTER_F_MODEL_H


#include "IExporter.h"
#include "../FrundFacade/FModel.h"


/**
* Класс, предназначенный для экспорта модели FRUND в различные форматы хранения данных
*/
class ExporterFModel
	:
		public IExporter
{
public:

	// Конструкторы и деструктор

	/**
	* Конструктор
	* @param model
	*/
	ExporterFModel
		(
			const FModel* model
		);


	// Селекторы

	/**
	* Получить наименование файла с расширением, в который будет производиться экспорт модели
	* @return наименование файла с расширением, в который будет производиться экспорт модели
	*/
	string GetModelFileName() const;


	// Функции экспорта

#pragma region overriden

	/**
	* Экспорт по умолчанию
	*/
	virtual
	void Export() const override;

#pragma endregion

	/**
	* Экспорт модели в указанную директорию
	* @param targetDirectory - директория для сохранения модели
	*/
	void ExportToDirectory
		(
			const string& directory
		)	const;

	/**
	* Экспорт модели в указанную директорию
	* @param targetDirectory - директория для сохранения модели
	* @param targetFile - файл для сохранения модели
	*/
	void ExportToDirectory(const string& targetDirectory, const string& targetFile) const;

protected:

	// ссылка на модель
	const FModel* _model;

};


#endif // EXPORTER_F_MODEL_H