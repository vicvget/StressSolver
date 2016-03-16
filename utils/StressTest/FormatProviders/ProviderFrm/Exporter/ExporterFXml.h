#ifndef EXPORTER_FXML_H

#define EXPORTER_FXML_H


#include "ExporterFModel.h"
#include "../../ProviderXml/XMLHelper.h"


/** Класс для экспорта FModel в XML формат
*
* @author Getmanskiy Victor
*/
class ExporterFXml :
	public ExporterFModel
{
public:

	ExporterFXml
		(
			const FModel* model
		);

#pragma region overriden

	/** Экспорт в заданный файл
	* @param file - файл
	*/
	virtual void ExportTo(const string& file) const override;

	virtual string GetExt() const override;

#pragma endregion


private:

	TiXmlElement* ExportHeader() const;

	TiXmlElement* ExportBody(size_t bodyId) const;
	TiXmlElement* ExportGeometry(size_t geometryId) const;
	TiXmlElement* ExportJoint(size_t id) const;
	TiXmlElement* ExportJointDofDescription(size_t dofId, CharComponent spring, CharComponent damping) const;
	TiXmlElement* ExportJointChar(size_t id) const;

	TiXmlElement* ExportForce(size_t id) const;
	TiXmlElement* ExportForceDofDescription(size_t dofId, CharComponent force) const;
	TiXmlElement* ExportForceChar(size_t id) const;

//	TiXmlElement* ExportForce(size_t id) const;
//	TiXmlElement* ExportForceChar(size_t id) const;

};


#endif // EXPORTER_FXML_H