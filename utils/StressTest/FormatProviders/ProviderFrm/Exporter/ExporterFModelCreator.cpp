#include "ExporterFModelCreator.h"

#include "ExporterFrm.h"
#include "ExporterFXml.h"
#include "ExporterNxXml.h"
#include "ExporterGraphViz.h"


/**
* ������� ������-��������� ���������� ������ ����� � ������ ������� ����
* @param model - ����������� ������ �����
* @param format - ������ ���������� ������ �����
* @return ��������� ������-��������� ���������� ������ �����
*/
ExporterFModel* CreateExporterFModel
	(
		const FModel* model,
		FModelFormat format
	)
{
	ExporterFModel* exporter = nullptr;

	switch (format)
	{
	case FModelFormat::frm:
		exporter = new ExporterFrm(model);
		break;

	case FModelFormat::xml:
		exporter = new ExporterFXml(model);
		break;

	case FModelFormat::nx:
		exporter = new ExporterNxXml(model);

	case FModelFormat::gv:
		exporter = new ExporterGraphViz(model);

	}

	return exporter;
}