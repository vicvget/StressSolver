#ifndef EXPORTER_F_MODEL_CREATOR_H

#define EXPORTER_F_MODEL_CREATOR_H


#include "ExporterFModel.h"


// ������ �������� ���������� ������
enum class FModelFormat
{
	frm,	// ���������� � ����� *.frm
	xml,	// ���������� � *.xml (������ FRUND)
	nx,		// ���������� � *.xml (������ Siemens NX)
	gv		// ���������� � *.gv (GraphViz)
};


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
	);


#endif // EXPORTER_F_MODEL_CREATOR_H