#ifndef EXPORTER_NX_XML_H

#define EXPORTER_NX_XML_H


#include "ExporterFModel.h"
#include "../FrundFacade/FModel.h"
#include "../../ProviderXml/XMLHelper.h"
#include "MechanismOccurrenceData.h"
#include <map>

using namespace MathHelpers;

/** Класс для экспорта FModel в XML формат Siemens NX
*
* @author Getmanskiy Victor
*/
class ExporterNxXml :
	public ExporterFModel
{
	MechanismOccurrence::MechanismOccurrenceData _mechanismOccurrenceData;
public:
	
	ExporterNxXml
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
	// добавляет элемент с id

	TiXmlElement* AddUserValue(const char* title, int value, TiXmlElement* parent) const;
	TiXmlElement* AddUserValue(const char* title, const char* value, TiXmlElement *parent) const;

	TiXmlElement* AddElement(int& id, const char* element, TiXmlElement* parent) const;
	TiXmlElement* AddElement(int& id, const char* element, const char* type, TiXmlElement* parent) const;


	TiXmlElement* AddMechanismRevisionView(int& id, const string& name, const string subtype, TiXmlElement* parent) const;
	TiXmlElement* AddMarker(int& id, const string& name, TiXmlElement* parent) const;
	TiXmlElement* AddMarker(int& id, const string& name, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent) const;
	
	TiXmlElement* AddTransform(int& id, const Vec3& translate, TiXmlElement* parent) const;
	TiXmlElement* AddTransform(int& id, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent) const;
	TiXmlElement* AddInertia(int& id, double mass, const Vec3& inertia, const Vec3& cmPoint, const Vec3& translate, const Mat3& rotation, TiXmlElement* parent, int& cmNodeId) const;

	TiXmlElement* ExportBody(int& id, int bodyId, TiXmlElement* parent) const;
	TiXmlElement* ExportJoint(int& id, int instancedRef, int jointId, TiXmlElement* parent) const;
	TiXmlElement* ExportForce(int& id, int forced, TiXmlElement* parent) const;

	TiXmlElement* ExportJointForces(int& id, int instancedRef, int jointId, TiXmlElement* parent) const;
	
	void ExportDampers(int& id, int instancedRef, TiXmlElement* parent) const;
	void ExportSprings(int& id, int instancedRef, TiXmlElement* parent) const;

	TiXmlElement* AddConstraintSeparator(int& id, const char* type, TiXmlElement* parent) const;

	TiXmlElement* ExportForces(int& id, int instancedRef, int forceId, TiXmlElement* parent, TiXmlElement* header) const;
	TiXmlElement* ExportTorques(int& id, int instancedRef, int forceId, TiXmlElement* parent, TiXmlElement* header) const;

	TiXmlElement* AddConstraintTargetRef(int& id, TiXmlElement* xmlConstraintInstance, int bodyId, int nodeId) const;
	TiXmlElement* AddConstraintTargetRefToMarker(int& id, TiXmlElement* xmlConstraintInstance, int markerId) const;


	struct BodyInfo
	{
		bool _isGround;
		int _cmNodeXmlId;
		std::vector<int> nodeIds;
	};

	struct MarkerInfo
	{
		Vec3 _translation;
		Mat3 _rotation;
		TiXmlElement* _xmlTransform;

		MarkerInfo()
			:
				_xmlTransform(nullptr)
		{
		}

		MarkerInfo(const Vec3& translate, const Mat3& rotation, TiXmlElement* xmlTransform):
			_translation(translate),
			_rotation(rotation),
			_xmlTransform(xmlTransform){}

		void AddTransformToXML() const;
	};


	struct SpringDamperInfo
	{
		std::string _name;
		double _stiffnessFactor;
		double _dampingFactor;
		double _preload;
		double _length;
		std::vector<int> _markerIds;
	};

	mutable std::map<int, MarkerInfo> _markerInfos;
	mutable std::vector<BodyInfo> _bodyInfos;
	mutable std::vector<SpringDamperInfo> _springDamperInfos;


};


#endif // EXPORTER_NX_XML_H