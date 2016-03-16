#include "ExporterFModel.h"


// class ExporterFModel

ExporterFModel::ExporterFModel
	(
		const FModel* model
	)
	:
		_model(model)
{
}

/**
* Получить наименование файла с расширением, в который будет производиться экспорт модели
* @return наименование файла с расширением, в который будет производиться экспорт модели
*/
string ExporterFModel::GetModelFileName() const
{
	return _model->Name() + GetExt();
}

/**
* Экспорт по умолчанию
*/
// virtual
void ExporterFModel::Export() const // override
{
	ExportTo(GetModelFileName());
}

/**
* Экспорт модели в указанную директорию
* @param targetDirectory - директория для сохранения модели
*/
// virtual
void ExporterFModel::ExportToDirectory
	(
		const string& targetDirectory
	)	const
{
	string targetPath = fs::CombinePath(targetDirectory, GetModelFileName());

	ExportTo(targetPath);
}

/**
* Экспорт модели в указанную директорию
* @param targetDirectory - директория для сохранения модели
*/
// virtual
void ExporterFModel::ExportToDirectory
(
const string& targetDirectory,
const string& targetFile
)	const
{
	string targetPath = fs::CombinePath(targetDirectory, targetFile);

	ExportTo(targetPath);
}