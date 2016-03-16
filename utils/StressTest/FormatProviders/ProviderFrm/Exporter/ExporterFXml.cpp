#include "ExporterFXml.h"

#include "../../../Fcore/fcore.h"


ExporterFXml::ExporterFXml
	(
		const FModel* model
	)
	:
		ExporterFModel(model)
{
}

TiXmlElement* ExporterFXml::ExportHeader() const
{
	const FHeader& header = _model->Header();
	XMLHelper xml;
	TiXmlElement* xmlBody = xml.InsertNewElement("header");
	xml.InsertConstElement("veriableGroupId", header.GetVariablesGroupsSet().CurrentId());
	for (unsigned int i = 0; i < header.GetVariablesGroupsCount(); i++)
	{
		xml.InsertNewElement("variableGroup", i);
		const FVariablesGroup& variableGroup = header.GetVariablesGroup(i);
		xml.InsertConstElement("name", variableGroup.Name());
		
		for (int j = 0; j < variableGroup.Size(); j++)
		{
			xml.InsertNewElement("variable");
			xml.InsertConstElement("name", variableGroup[j].name);
			xml.InsertConstElement("expression", variableGroup[j].expression);
			xml.InsertConstElement("description", variableGroup[j].description);
			xml.LevelUp();
		}

		xml.LevelUp();
	}

	// TODO: варианты начальных условий, параметров и т.д.
	return xmlBody;
}


TiXmlElement* ExporterFXml::ExportBody(size_t bodyId) const
{
	const FBody& body = _model->GetBody(bodyId);
	XMLHelper xml;
	TiXmlElement* xmlBody = xml.InsertNewElement("body",body.Number());
	xml.InsertConstElement("name", body.Name());
	xml.InsertConstElement("number", body.Number());
	xml.InsertNewElement("inertia");
	xml.InsertNewElement("standart");
	xml.InsertConstElement("mass",body.M());
	xml.InsertNewElement("inertiaTensor");
	xml.InsertConstElement("Ixx",body.Jx());
	xml.InsertConstElement("Iyy",body.Jy());
	xml.InsertConstElement("Izz",body.Jz());
	xmlBody->InsertEndChild(*ExportGeometry(bodyId));
	return xmlBody;
}


TiXmlElement* ExporterFXml::ExportForce(size_t id) const
{
	const FForce& force = _model->GetForce(id);
	XMLHelper xml;
	TiXmlElement* xmlForce = xml.InsertNewElement("force", force.Number());
	xml.InsertConstElement("name", force.Name());
	xml.InsertConstElement("number", force.Number());
	xml.InsertRefElement("forceChar", force.CharNumber());
	xml.InsertNewElement("bodyPoint", 0);
	xml.InsertRefElement("body", force.BodyNumber());
	xml.InsertRefElement("point", force.NodeNumber());
	return xmlForce;
}

TiXmlElement* ExporterFXml::ExportForceDofDescription(size_t dofId, CharComponent forceComponent) const
{
	XMLHelper xml;
	TiXmlElement* xmlDof = xml.InsertNewElement("dof", dofId);

	XMLHelper xmlSpring;
	xml.InsertConstElement("type", forceComponent.type);
	xml.InsertConstElement("direction", forceComponent.direction);
	xml.InsertNewElement("params");
	for (int i = 0; i < forceComponent.params.size(); i++)
	{
		xml.InsertConstElement("param", forceComponent.params[i]);
	}
	return xmlDof;
}

TiXmlElement* ExporterFXml::ExportForceChar(size_t id) const
{
	const FForceChar& forceChar = _model->GetForceChar(id);
	XMLHelper xml;
	TiXmlElement* xmlForceChar = xml.InsertNewElement("forceChar", forceChar.Number());
	xml.InsertConstElement("name", forceChar.Name());
	xml.InsertConstElement("number", forceChar.Number());

	auto forceDofs = forceChar.GetParams();
	const int nDofs = forceDofs.size();

	for (int i = 0; i < nDofs; i++)
	{
		xmlForceChar->InsertEndChild(*ExportForceDofDescription(i, forceDofs[i]));
	}

	return xmlForceChar;
}


TiXmlElement* ExporterFXml::ExportJointDofDescription(size_t dofId, CharComponent springComponent, CharComponent dampingComponent) const
{
	XMLHelper xml;
	TiXmlElement* xmlDof = xml.InsertNewElement("dof", dofId);

	XMLHelper xmlSpring;
	xml.InsertNewElement("spring", dofId);
	xml.InsertConstElement("type", springComponent.type);
	xml.InsertConstElement("direction", springComponent.direction);
	xml.InsertNewElement("params");
	for (int i = 0; i < springComponent.params.size(); i++)
	{
		xml.InsertConstElement("param", springComponent.params[i]);
	}
	xml.LevelUp();
	xml.LevelUp();

	TiXmlElement* dampingNode = xml.InsertNewElement("damping", dofId);

	// TODO: использовать переменную dampingNode
	(void)dampingNode;
	xml.InsertConstElement("type", dampingComponent.type);
	xml.InsertConstElement("direction", dampingComponent.direction);
	xml.InsertNewElement("params");
	for (int i = 0; i < dampingComponent.params.size(); i++)
	{
		xml.InsertConstElement("param", dampingComponent.params[i]);
	}
	return xmlDof;
}


TiXmlElement* ExporterFXml::ExportJointChar(size_t jointCharId) const
{
	const FJointChar& jointChar = _model->GetJointChar(jointCharId);
	XMLHelper xml;
	TiXmlElement* xmlJointChar = xml.InsertNewElement("jointChar", jointChar.Number());
	xml.InsertConstElement("name", jointChar.Name());
	xml.InsertConstElement("number", jointChar.Number());

	auto springDofs = jointChar.GetSpringParams();
	auto dampingDofs = jointChar.GetDampingParams();
	const int nDofs = springDofs.size();

	for (int i = 0; i < nDofs; i++)
	{
		xmlJointChar->InsertEndChild(*ExportJointDofDescription(i, springDofs[i], dampingDofs[i]));
	}

	return xmlJointChar;
}


TiXmlElement* ExporterFXml::ExportJoint(size_t jointId) const
{
	const FJoint& joint = _model->GetJoint(jointId);
	XMLHelper xml;
	TiXmlElement* xmlJoint = xml.InsertNewElement("joint", joint.Number());
	xml.InsertConstElement("name", joint.Name());
	xml.InsertConstElement("number", joint.Number());
	xml.InsertRefElement("jointChar", joint.CharNumber());
	xml.InsertNewElement("bodyPoint", 0);
	xml.InsertRefElement("body", joint.BodyNodeNumber1().bodyNumber);
	xml.InsertRefElement("point", joint.BodyNodeNumber1().nodeNumber);
	xml.LevelUp();
	xml.InsertNewElement("bodyPoint", 1);
	xml.InsertRefElement("body", joint.BodyNodeNumber1().bodyNumber);
	xml.InsertRefElement("point", joint.BodyNodeNumber1().nodeNumber);
	return xmlJoint;
}


TiXmlElement* ExporterFXml::ExportGeometry(size_t geometryId) const
{
	const FGeometry& geo = _model->GetBodyGeometry(geometryId);
	XMLHelper xml;
	TiXmlElement* xmlGeo = xml.InsertNewElement("geometry");
	xml.InsertRefElement("centerOfMass",geo.CmNodeNumber());
	for(int i = 0; i < geo.CountNodes(); i++)
	{
		FGeometryPoint node = geo.GetNode(i);
		xml.InsertNewElement("point",i);
		xml.InsertConstElement("x",node.x);
		xml.InsertConstElement("y",node.y);
		xml.InsertConstElement("z",node.z);
		xml.LevelUp();
	}

	for(int i = 0; i < geo.CountLinks(); i++)
	{
		FGeometryLink link = geo.GetLink(i);
		xml.InsertNewElement("link",i);
		xml.InsertRefElement("point", link.node1);
		xml.InsertRefElement("point", link.node2);
		xml.LevelUp();
	}
	return xmlGeo;
}

void ExporterFXml::ExportTo(const string& file) const
{
	TiXmlDocument doc;  
	doc.LinkEndChild(new TiXmlDeclaration("1.0", "", ""));  
	TiXmlElement * root = new TiXmlElement("task");  
	doc.LinkEndChild(root);  

	root->LinkEndChild(new TiXmlComment("Model Description"));  
	TiXmlElement *xmlModel = new TiXmlElement("model");
	root->LinkEndChild(xmlModel);  
	xmlModel->SetAttribute("id", _model->Header().Name().c_str());
	// TODO: экспорт заголовка


	xmlModel->LinkEndChild(ExportHeader());

	// экспорт тел
	for(int i = 0; i<_model->CountBodies(); i++)
	{
		xmlModel->LinkEndChild(ExportBody(i));
	}
	
	// экспорт характеристик связей
	for (int i = 0; i<_model->GetJointCharsCount(); i++)
	{
		std::cout << "Char " << i << " exported\n";
		xmlModel->LinkEndChild(ExportJointChar(i));
	}
	
	// экспорт связей
	for (int i = 0; i<_model->GetJointsCount(); i++)
	{
		xmlModel->LinkEndChild(ExportJoint(i));
	}
	
	// экспорт характеристик сил
	for (int i = 0; i<_model->GetForceCharsCount(); i++)
	{
		xmlModel->LinkEndChild(ExportForceChar(i));
	}
	
	// экспорт сил
	for (int i = 0; i<_model->GetForcesCount(); i++)
	{
		xmlModel->LinkEndChild(ExportForce(i));
	}

	doc.SaveFile(file.c_str());
}

string ExporterFXml::GetExt() const
{
	return EXT_XML;
}
