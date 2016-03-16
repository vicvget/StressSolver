#ifndef EXPORTER_F_MODEL_H

#define EXPORTER_F_MODEL_H


#include "IExporter.h"
#include "../FrundFacade/FModel.h"


/**
* �����, ��������������� ��� �������� ������ FRUND � ��������� ������� �������� ������
*/
class ExporterFModel
	:
		public IExporter
{
public:

	// ������������ � ����������

	/**
	* �����������
	* @param model
	*/
	ExporterFModel
		(
			const FModel* model
		);


	// ���������

	/**
	* �������� ������������ ����� � �����������, � ������� ����� ������������� ������� ������
	* @return ������������ ����� � �����������, � ������� ����� ������������� ������� ������
	*/
	string GetModelFileName() const;


	// ������� ��������

#pragma region overriden

	/**
	* ������� �� ���������
	*/
	virtual
	void Export() const override;

#pragma endregion

	/**
	* ������� ������ � ��������� ����������
	* @param targetDirectory - ���������� ��� ���������� ������
	*/
	void ExportToDirectory
		(
			const string& directory
		)	const;

	/**
	* ������� ������ � ��������� ����������
	* @param targetDirectory - ���������� ��� ���������� ������
	* @param targetFile - ���� ��� ���������� ������
	*/
	void ExportToDirectory(const string& targetDirectory, const string& targetFile) const;

protected:

	// ������ �� ������
	const FModel* _model;

};


#endif // EXPORTER_F_MODEL_H