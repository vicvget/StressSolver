#pragma once
#include "ExporterFModel.h"


/**  ласс дл€ экспорта FModel в GV формат
*
* @author Getmanskiy Victor
*/
class ExporterGraphViz: public ExporterFModel
{
public:

	ExporterGraphViz(const FModel* model);

#pragma region overriden

	/** Ёкспорт в заданный файл
	* @param file - файл
	*/
	virtual void ExportTo(const string& file) const override;

	virtual string GetExt() const override;

#pragma endregion


private:
};