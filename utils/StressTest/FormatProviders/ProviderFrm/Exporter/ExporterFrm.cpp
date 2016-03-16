#include "ExporterFrm.h"

#include "../../../Fcore/fcore.h"


ExporterFrm::ExporterFrm
	(
		const FModel* model
	)
	:
		ExporterFModel(model)
{
}

void ExporterFrm::ExportTo(const string& file) const
{
	_model->SaveTo(file);
}

string ExporterFrm::GetExt() const
{
	return EXT_FRM;
}