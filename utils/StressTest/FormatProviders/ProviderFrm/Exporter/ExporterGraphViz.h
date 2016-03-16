#pragma once
#include "ExporterFModel.h"


/** ����� ��� �������� FModel � GV ������
*
* @author Getmanskiy Victor
*/
class ExporterGraphViz: public ExporterFModel
{
public:

	ExporterGraphViz(const FModel* model);

#pragma region overriden

	/** ������� � �������� ����
	* @param file - ����
	*/
	virtual void ExportTo(const string& file) const override;

	virtual string GetExt() const override;

#pragma endregion


private:
};