#ifndef EXPORTER_F_MODEL_CREATOR_H

#define EXPORTER_F_MODEL_CREATOR_H


#include "ExporterFModel.h"


// Список форматов сохранения модели
enum class FModelFormat
{
	frm,	// сохранение в файлы *.frm
	xml,	// сохранение в *.xml (формат FRUND)
	nx,		// сохранение в *.xml (формат Siemens NX)
	gv		// сохранение в *.gv (GraphViz)
};


/**
* Создать объект-экспортер сохранения модели ФРУНД в формат данного типа
* @param model - сохраняемая модель ФРУНД
* @param format - формат сохранения модели ФРУНД
* @return созданный объект-экспортер сохранения модели ФРУНД
*/
ExporterFModel* CreateExporterFModel
	(
		const FModel* model,
		FModelFormat format
	);


#endif // EXPORTER_F_MODEL_CREATOR_H