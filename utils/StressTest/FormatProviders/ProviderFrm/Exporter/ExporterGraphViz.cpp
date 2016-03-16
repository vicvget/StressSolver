#include "ExporterGraphViz.h"
#include <fstream>

#include "../../../Fcore/fcore.h"

using std::ofstream;

ExporterGraphViz::ExporterGraphViz(const FModel* model) : ExporterFModel(model)
{

}

string ExporterGraphViz::GetExt() const
{
	return EXT_GV;
}

/** Ёкспорт в заданный файл
* @param file - файл
*/
void ExporterGraphViz::ExportTo(const string& file) const
{
	ofstream ofs(file);

	const string& stiffJointAttribute = "[color = red, penwidth = 5.0]";

	if (ofs.is_open())
	{
		ofs << "graph {\n";
		for (int i = 0; i < _model->GetJointsCount(); i++)
		{
			const FJoint& joint = _model->GetJoint(i);
			const FJointChar& jointChar = _model->GetJointCharByJoint(joint);


			bool isStiffJoint = false;
			if (jointChar.GetSpringParams().front().type == 15)
				isStiffJoint = true;
			const string& attribute = isStiffJoint ? stiffJointAttribute : string();

			ofs << joint.BodyNumber1() << " -- " << joint.BodyNumber2() << attribute << ";" << std::endl;
		}


		ofs << "}\n";
		ofs.close();
	}
}
