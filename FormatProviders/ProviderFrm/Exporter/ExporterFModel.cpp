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
* �������� ������������ ����� � �����������, � ������� ����� ������������� ������� ������
* @return ������������ ����� � �����������, � ������� ����� ������������� ������� ������
*/
string ExporterFModel::GetModelFileName() const
{
	return _model->Name() + GetExt();
}

/**
* ������� �� ���������
*/
// virtual
void ExporterFModel::Export() const // override
{
	ExportTo(GetModelFileName());
}

/**
* ������� ������ � ��������� ����������
* @param targetDirectory - ���������� ��� ���������� ������
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
* ������� ������ � ��������� ����������
* @param targetDirectory - ���������� ��� ���������� ������
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