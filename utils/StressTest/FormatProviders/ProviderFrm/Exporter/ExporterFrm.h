#ifndef EXPORTER_FRM_H

#define EXPORTER_FRM_H


#include "ExporterFModel.h"
#include "../FrundFacade/FModel.h"


class ExporterFrm :
	public ExporterFModel
{
public:

	ExporterFrm
		(
			const FModel* model
		);

#pragma region overriden

	virtual void ExportTo(const string& file) const override;

	virtual string GetExt() const override;

#pragma endregion

};


#endif // EXPORTER_FRM_H