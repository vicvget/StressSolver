#include "FModel.h"

#include "../../../fcore/fcore.h"
#include "../../../Fcore/wrappers/DirectoryLister.h"
#include "FMacroSpecification.h"

#include <algorithm>
#include <string>
#include <list>
#include <map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>


using Calc::Calculator;
using std::list;
using std::ios;
using std::setw;


char blockHeaders[][256] = {
	"BL31 NPDS    MASSA       JX       JY       JZ\n",			//0
	"BL14 NPDS NTIP   LX   LY   LZ   UX   UY   UZ\n",			//1
	"BL33  NMXAR    PNAPR    PSLAG     KPAR\n",					//2
	"BL16  NMX NAPR  TUF KPAR  TDF KPAR\n",						//3
	"BL17  NSE  NT1 NUZ1  NT2 NUZ2  NMX\n",						//4
	"BL34  NSS     NAPR     KSTR\n",							//5
	"BL18  NSS NAPR  NTF KPAR\n",								//6
	"BL19  NSS   NT  NUZ\n",									//7
	"BL32 NPDS      NUZ        X        Y        Z\n",			//8
	"BL35 NPDS   KSTROK     KPAR      ID1      ID2 \n",			//9
	"BL36 NPDS1   NUZL1    NPDS2    NUZL2     SKOL\n",			//10
	"BL37 NKOL       HS        D       FZ       D1       D2\n"	//11
};

// TODO:
/*
9 4 4 4 ...
"BL14 NPDS NTIP   LX   LY   LZ   UX   UY   UZ\n",
"BL16  NMX NAPR  TUF KPAR  TDF KPAR\n",
"BL17  NSE  NT1 NUZ1  NT2 NUZ2  NMX\n",
"BL18  NSS NAPR  NTF KPAR\n",
"BL19  NSS   NT  NUZ\n",


9 8 8 8 ...
"BL31 NPDS    MASSA       JX       JY       JZ\n",
"BL32 NPDS      NUZ        X        Y        Z\n",
"BL33  NMX    PNAPR    PSLAG     KPAR\n",
"BL34  NSS     NAPR     KSTR\n",
"BL35 NPDS   KSTROK     KPAR      ID1      ID2 \n",
"BL36 NPDS     NUZL    NPDS2    NUZL2     SKOL\n",
"BL37 NKOL       HS        D       FZ       D1       D2\n"
*/

/**
* Формирование карты конверсии номеров тел в индексы в массиве bodies
*/
void FModel::FormNumbersToIndexes()
{
	_numbersToIndexes.clear();
	for (int id = 0; id < GetBodiesCount(); id++)
	{
		_numbersToIndexes[_bodies[id].Number()] = id;
	}
}

void FModel::Save() const
{
	Save(GetModelFileName());
}

void FModel::SaveTo(const string& newPath) const
{
	if(newPath.empty())
		Save();
	else
	{
		string dir = fs::SplitDirFromPath(newPath);
		string file = fs::SplitFileFromPath(newPath);

		if(file.empty())
			file = _name;
		
		if(dir.empty())
		{
			Save(file);
		}
		else
		{
			string oldPath;

			fs::MakeDir(dir);
			// debug 
			//std::cout << fs::GetCurrentDir() << " after make" << std::endl;
			vector<string> filesToCopy;
			filesToCopy.push_back("*.ico");
			filesToCopy.push_back("*.bnf");
			filesToCopy.push_back("*.flag");
			filesToCopy.push_back("*.rlc");
			filesToCopy.push_back("*.stl");
			filesToCopy.push_back("*.stp");
			filesToCopy.push_back("*.step");
			filesToCopy.push_back("fadres.dat");
			fs::CopyFilesToFolder(filesToCopy, dir);
			fs::SetCurrentDir(dir,oldPath);
			Save(file);
			fs::SetCurrentDir(oldPath);
		}
	}
}

void FModel::Save(const string& path) const
{
	_debugParams.Save();
	ofstream stream(path);
	if(!stream.is_open())
		exceptions::ThrowFileNotOpened(path);
	Header().SaveByIndex();
	stream << 3 << ' ' << Header().GetFileIndex() << std::endl;
	size_t i;
	for(i = 0; i < GetBodiesCount(); i++)
	{
		_bodies[i].SaveByIndex();
		stream << 1 << ' ' << _bodies[i].GetFileIndex() << std::endl;
	}
	for(i = 0; i < GetGeometriesCount(); i++)
	{
		_geometries[i].SaveByIndex();
	}
	for(i = 0; i < GetJointCharsCount(); i++)
	{
		_jointChars[i].SaveByIndex();
		stream << 2 << ' ' << _jointChars[i].GetFileIndex() << std::endl;
	}
	for(i = 0; i < GetCunitsCount(); i++)
	{
		_cunits[i].SaveByIndex();
		stream << 4 << ' ' << _cunits[i].GetFileIndex() << std::endl;
	}
	for(i = 0; i < GetJointsCount(); i++)
	{
		_joints[i].SaveByIndex();
		stream << 5 << ' ' << _joints[i].GetFileIndex() << std::endl;
	}
	for(i = 0; i < GetForceCharsCount(); i++)
	{
		_forceChars[i].SaveByIndex();
		stream << 6 << ' ' << _forceChars[i].GetFileIndex() << std::endl;
	}
	for(i = 0; i < GetForcesCount(); i++)
	{
		_forces[i].SaveByIndex();
		stream << 7 << ' ' << _forces[i].GetFileIndex() << std::endl;
	}
	for (i = 0; i < GetFreeSolversCount(); i++)
	{
		_freeSolvers[i].SaveByIndex();
		stream << 11 << ' ' << _freeSolvers[i].GetFileIndex() << std::endl;
	}
	stream.close();


	Header().SaveByIndex();
	
	stream.close();

}

FModel::FModel(const char* path)
{
	Load(path);
}

FModel::FModel() :_name("default.frm")
{

}

void FModel::Load(const char* path)
{
	string index;

	if (!fs::FileExist(path))
	{
		exceptions::ThrowFileNotFound(path);
	}
	_name = fs::SplitFileFromPath(path);
	_name = fs::SplitFileFromExt(_name);

	ifstream stream(path);

	if(!stream.is_open()) 
		exceptions::ThrowFileNotOpened(path);
	else
	{
		_debugParams.Load();

		//_path = path;

		int type = 0;

		stream >> type >> index;
		while(!stream.eof())
		{
			index = index.append(EXT_MODEL_ELEMENT);
			switch(type)
			{
			case 1: // body
				{
					FBody body(index.c_str());
					body.Id(GetBodiesCount());
					FGeometry geometry;
					geometry.Init(body.FileId(),EXT_BODY_GEOMETRY,body.Id(), body.Number());
					if(geometry.LoadByIndex())
					{
						body.GeometryId(GetGeometriesCount()); // 0-based ?
						_geometries.push_back(geometry);
					}
					_bodies.push_back(body);
				}
				break;
			case 2: // joint char
				{
					FJointChar jointChar(index.c_str());
					_jointChars.push_back(jointChar);
				}
				break;
			case 3: // header
				{
					//header.FileId(index.c_str());
					_header.Load(index.c_str());
				}
				break;
			case 4: // CU ???
				{
					FCunit cunit(index.c_str());
					_cunits.push_back(cunit );
				}
				break;
			case 5: // joint
				{
					FJoint joint(index.c_str());
					_joints.push_back(joint);
				}
				break;
			case 6: // force char
				{
					FForceChar forceChar(index.c_str());
					_forceChars.push_back(forceChar);
				}
				break;
			case 7: // force
				{
					FForce force(index.c_str());
					_forces.push_back(force);
				}
				break;
			case 11: // free solver
				{
					string new_index = index.substr(0, index.find_last_not_of(EXT_MODEL_ELEMENT) + 1);
					FFreeSolver freeSolver(new_index.c_str());
					_freeSolvers.push_back(freeSolver);
				}
				break;

			default:
				break;
			}
			stream >> type >> index;
		}
		stream.close();
		FormNumbersToIndexes();
	}
}

/**
* Загрузить геометрию профиля из *949.gmr или *49.gmr файлов
* @return признак успешной (true) или неуспешной (false) загрузки
*/
bool FModel::LoadRoadProfileGeometry()
{
	if (!IsBodyExistByNumber(49))
	{
		return false;
	}

	string fileRoadGmr;
	FBody roadBody = GetBodyByNumber(49);
	ifstream ifsRoad;

	fileRoadGmr = Header().GetDetiledRoadProfile();
	if(!fs::FileExist(fileRoadGmr))
		fileRoadGmr = Header().GetRoadProfile();
	if(!fs::FileExist(fileRoadGmr))
		return false;

	ifsRoad.open(fileRoadGmr.c_str());
	if (ifsRoad.is_open())
	{
		FGeometry roadGeometry(ifsRoad, 1000.);
		if(roadBody.HasGeometry())
		{
			_geometries[roadBody.GeometryId()] = roadGeometry;
		}
		else
		{
			roadBody.GeometryId(GetGeometriesCount());
			_geometries.push_back(roadGeometry);
		}
		_bodies[roadBody.Id()] = roadBody;
		ifsRoad.close();
		return true;
	}
	return false;
}

/**
* Получить наименование файла модели с расширением
* @return наименование файла модели с расширением
*/
string FModel::GetModelFileName() const
{
	return _name + EXT_FRM;
}

/**
* Возвращает индекс тела в массиве bodies по его номеру
* @param bodyNumber - номер тела
* @return индекс тела в массиве bodies (-1, если тело с таким номером не найдено)
*/
int FModel::GetBodyIndexByNumber
	(
		int bodyNumber
	)	const
{
	IndexConversionMap::const_iterator it;

	it = _numbersToIndexes.find(bodyNumber); 
	if (it != _numbersToIndexes.end())
	{
		return it->second;
	}
	else
	{
		// TODO: исправить на выброс исключения

		return -1; 
	}
}

// Вывод расчетной схемы в формате ФРУНД
void FModel::Output(ofstream& stream, int isRType)
{
	Eval();
	OutputHeader(stream, isRType);
	OutputVariables(stream, isRType);
	OutputBodies(stream, isRType);
	OutputCunits(stream, isRType);
	OutputJointChars(stream, isRType);
	if(!isRType) 
		SortOutputSequence(3);
	SortOutputSequence(4);
	OutputJoints(stream, isRType);
	OutputForceChars(stream, isRType);
	if(!isRType) OutputForces(stream, isRType);
	if(isRType) OutputAdditionalVehicleParams(stream);
	OutputMacros(stream, isRType);
	stream << "ENDE\n";
}

bool comp(const FJoint& el1, const FJoint& el2)
{
	return el1.SortId() < el2.SortId();
}

class FindGearCharPred
{
public:

	FindGearCharPred
		(
			const FModel& model,
			const FJoint& joint
		)
		:
			_model(model),
			_joint(joint)
	{
	}

	bool operator()(const FJoint& joint2)
	{
		if (_joint.Number() == joint2.Number())
			return false;

		int type1, type2;
		size_t subtype1, subtype2;

		_model.GetJointCharByJoint(_joint).GetSpringParamsType(0, type1, subtype1);
		_model.GetJointCharByJoint(joint2).GetSpringParamsType(0, type2, subtype2);

		if (type1 == type2 && subtype1 == subtype2)
		{
			if
				(
					_joint.BodyNodeNumber1() == joint2.BodyNodeNumber1() ||
					_joint.BodyNodeNumber1() == joint2.BodyNodeNumber2() ||
					_joint.BodyNodeNumber2() == joint2.BodyNodeNumber1() ||
					_joint.BodyNodeNumber2() == joint2.BodyNodeNumber2()
				)
			{
				return true;
			}
		}

		return false;
	}

private:

	const FModel& _model;

	const FJoint& _joint;

};

// Сформировать последовательность вывода в sequence
// type = 0 - для тел, sequence = bodySequence[rest]
// type = 1 - для тел, по кол-ву уравнений связей (bl[0]=-1)
// type = 2 - для тел, по кол-ву уравнений связей (bl[0]=-2)
// type = 3 - для с.е., в обратном порядке следования тел bodySequence
// type = 4 - для с.е. перенос в конец характеристик шестеренок
void FModel::SortOutputSequence(int type)
{
	switch (type)
	{
	case 0:
		// TODO:
		break;

	case 1:
		// TODO:
		break;

	case 2:
		// TODO:
		break;

	case 3:
		{
			vector<int> bodySequence;

			for (const FBody& body : _bodies)
			{
				bodySequence.push_back(body.Number());
			}

			list<const FJoint*> pjoints;

			for (const FJoint& joint : _joints)
			{
				joint.SortId(std::numeric_limits<int>::max());
				pjoints.push_back(&joint);
			}

			list<const FJoint*>::const_iterator itPJoint, itRemove;
			int sortId = 0;

			for (int bodyNumber : bodySequence)
			{
				itPJoint = pjoints.begin();
				while (itPJoint != pjoints.end())
				{
					if
						(
							((*itPJoint)->BodyNumber1() == bodyNumber) ||
							((*itPJoint)->BodyNumber2() == bodyNumber)
						)
					{
						itRemove = itPJoint;
						(*itPJoint)->SortId(sortId++);
						++itPJoint;
						pjoints.remove(*itRemove);
					}
					else
					{
						++itPJoint;
					}
				}
			}
			sort
				(
					_joints.begin(),
					_joints.end(),
					[]
						(
							const FJoint& joint1,
							const FJoint& joint2
						)
					{
						return joint1.SortId() < joint2.SortId();
					}
				);
			reverse(_joints.begin(), _joints.end());
		}
		break;

	case 4:
		{
// Процедура поиска шестеренок и переноса их в конец в порядке (жёстк=0, жёстк=0.0002)
			Calculator* calc = GetCalculator();
			map<const FJoint*, const FJoint*> gearJoints;
			vector<FJoint>::const_iterator itJoint = _joints.begin();
			int type = 0;
			size_t subtype = 0;

			while (itJoint != _joints.end())
			{
				GetJointCharByJoint(*itJoint).GetSpringParamsType(0, type, subtype);
				if(type == 15 && subtype == 4)//gear
				{
					double stiffness = calc->Eval(GetJointCharByJoint(*itJoint).GetSpringParam(0,0).c_str());
					if(fabs(stiffness) < 1e-15)
					{
						// find adjacent joint
						FindGearCharPred pred(*this, *itJoint);
						vector<FJoint>::const_iterator itJoint2 = find_if(_joints.begin(),_joints.end(), pred);
						if(itJoint2 != _joints.end())
						{
							gearJoints[&(*itJoint)] = &(*itJoint2);
						}
						else
						{
							exceptions::ThrowNoElementInCollection("joint 15 4 pair for gears", "joints");
							// TODO: ошибка!
						}
					}
				}
				++itJoint;
			}
			delete calc;
			if (!gearJoints.empty())
			{
				vector<FJoint>::const_iterator itJoint = _joints.begin();
				int sortId = 0;
				while (itJoint != _joints.end())
				{
					itJoint->SortId(sortId++);
					++itJoint;
				}

				map<const FJoint*, const FJoint*>::const_iterator itGearJoint = gearJoints.begin();
				while(itGearJoint != gearJoints.end())
				{
					itGearJoint->first->SortId(sortId++);
					itGearJoint->second->SortId(sortId++);
					++itGearJoint;
				}

				sort(_joints.begin(), _joints.end(), comp);
			}
		}
		break;

	}
}


// Вывод в текстовый файловый поток заголовка расчетной схемы
void FModel::OutputHeader(ofstream& stream, int isRType) const
{
	_header.Output(stream, isRType);
}

// Вывод в текстовый файловый поток переменных текущего набора
void FModel::OutputVariables(ofstream &stream, int isRType) const
{
	_header.OutputVariables(stream, isRType);
}

// Вывод в текстовый файловый поток макросов
void FModel::OutputMacros(ofstream &stream, int isRType) const
{
	_header.OutputMacros(stream, isRType);
}

// Вывод информации о телах в модели
void FModel::OutputBodies(ofstream& stream, int isRType) const
{	
	int blockId = isRType ? 0 : 1;
	stream << blockHeaders[blockId];
	for(size_t i = 0; i < _bodies.size(); i++)
	{
		_bodies[i].MinMass(_debugParams.minMass);
		_bodies[i].MinInertia(_debugParams.minInertia);
		_bodies[i].Output(stream, isRType);
	}
}

// Вывод информации о фиктивных телах в модели
void FModel::OutputCunits(ofstream& stream, int isRType) const
{	
	for(size_t i = 0; i < _cunits.size(); i++)
	{
		_cunits[i].Output(stream, isRType);
	}
}

// Вывод информации о характеристиках шарниров в модели
void FModel::OutputJointChars(ofstream& stream, int isRType) const
{
	if(isRType) 
		fs::DeleteFiles(FILE_RHL_SPRINGS);
	else
		fs::DeleteFiles(FILE_MHL_SPRINGS);
	int blockId = isRType ? 2 : 3;
	stream << blockHeaders[blockId];
	for(size_t i = 0; i < _jointChars.size(); i++)
	{
		_jointChars[i].SetMaxStiffness(_debugParams.maxStiffness);
		_jointChars[i].Output(stream, isRType);
	}
}

// Вывод информации о шарнирах в модели
void FModel::OutputJoints(ofstream& stream, int isRType) const
{
	stream << blockHeaders[4];
	for(size_t i = 0; i < _joints.size(); i++)
	{
		_joints[i].Output(stream, isRType);
	}
}

// Вывод информации о характеристиках сил в модели
void FModel::OutputForceChars(ofstream& stream, int isRType) const
{
	int blockId = isRType ? 5 : 6;
	stream << blockHeaders[blockId];
	for(size_t i = 0; i < _forceChars.size(); i++)
	{
		_forceChars[i].Output(stream, isRType);
	}
}

//bool Equal(const FForce& force1, const FForce& force2)

/* function object to check the value of a map element
*/
//template <class K, class V>
class force_number_equals
{
public:

	// constructor (initialize value to compare with)
	force_number_equals(int v)
		: value(v)
	{
	}

	// comparison
	bool operator() (const FForce& force) const
	{
		return force.Number() == value;
	}

private:

	int value;

};

// Вывод информации о силах в модели
void FModel::OutputForces(ofstream& stream, int isRType) const
{
	stream << blockHeaders[7];
	for (int i = 0; i < _forceChars.size(); i++)
	{
		vector<FForce>::const_iterator it =
			find_if(_forces.begin(), _forces.end(), force_number_equals(_forceChars[i].Number()));

		if (it != _forces.end())
		{
			it->Output(stream, isRType);
		}
	}
}

// Вывод информации о телах в формате библиотеки ФРУНД
void FModel::SaveLba() const
{
	for(int i = 0; i < _bodies.size(); i++)
	{
		ofstream ofs(_bodies[i].GetPath(EXT_LBA));
		if(ofs.is_open())
		{
			_bodies[i].OutputLba(ofs);
			GetBodyGeometry(i).OutputLba(ofs);
			ofs.close();
		}
	}
}

//void FModel::SaveLbm(AvmodelInterface& avmodel, char* path)
//{
//	if(avmodel.InitLibrary(path))
//	{
//		for(int i = 0; i < bodies.size(); i++)
//		{
//			if(!bodies[i].SaveLbm(avmodel))
//			{
//				//TODO: exception
//				cout << "Body is not inserted\n";// << bodies[i].id << bodies[i].name << "не добавлено\n";
//			}
//		}		
//		avmodel.CloseLibrary();
//	}	
//}

void FModel::GetGraphCSR( int *xadj ,int *adjncy,int *vwgt )
{
	int d = 0;
	int x = 0;
	//перебираем все тела
	for (int bodyIndex=0; bodyIndex < _bodies.size(); bodyIndex++)
	{
		if(bodyIndex==0)
			xadj[bodyIndex] = 1;

		vwgt[bodyIndex] = 0;

		//перебираем все связи
		for (int jointIndex = 0; jointIndex<_joints.size(); jointIndex++)
		{
			if (_joints[jointIndex].BodyNumber1() == _bodies[bodyIndex].Number()) 
			{
				//adjncy[x] = bodyIndex;
				for (int i=0; i < _bodies.size();i++)
					if(_bodies[i].Number() == _joints[jointIndex].BodyNumber2())
						adjncy[x] = _bodies[i].Id();

				x++;
				d++;
				vwgt[bodyIndex]++;
			}
			else if (_joints[jointIndex].BodyNumber2() == _bodies[bodyIndex].Number())
			{
				//adjncy[x] = bodyIndex;

				for (int i=0; i < _bodies.size();i++)
					if(_bodies[i].Number() == _joints[jointIndex].BodyNumber1())
						adjncy[x] = _bodies[i].Id();
				x++;
				d++;
				vwgt[bodyIndex]++;
			}
		}
		xadj[bodyIndex+1] = xadj[bodyIndex]+d;
		d = 0;
	}
}

void FModel::OutputModelToGraph(ofstream& stream)
{
	stream << "Graph G { \n";
	for (size_t i = 0; i < _joints.size(); i++)
	{
		stream << _joints[i].BodyNumber1() << " -- " << _joints[i].BodyNumber2() << ";\n";
	}
	stream << "}\n";
}

// Вывод управляющих параметров для расчета 
// (uprf, epailon.dat, initcond.inp, way.cnt, tire.cnt)
void FModel::OutputSolveParams() const
{
	_header.OutputSolveParams(IsRoadExist());
}

// Вывод управляющих параметров для расчета (UPRF)
void FModel::OutputSolveControlParams(ofstream &stream) const
{
	_header.OutputSolveControlParams(stream);
}

// Вывод управляющих параметров для анализа (UPRF)
void FModel::OutputAnazControlParams(ofstream &stream) const
{
	_header.OutputAnazControlParams(stream);
}

/** Вывод начальных условий, соответствующих исходной геометрии
*/
void FModel::OutputDefaultIco() const
{
	ofstream ofs(FILE_DEFAULT_INITIAL_CONDITIONS);

	if (ofs.is_open())
	{
		for (int id = 0; id < GetBodiesCount(); id++)
		{
			int bodyNumber = _bodies[id].Number();
			ofs << bodyNumber << std::endl;

			const FGeometry& geo = GetBodyGeometry(id);
			ofs.precision(10);
			ofs.setf(ios::fixed);

			const float (&matrix)[12] = geo.GetTransformArray();

			if(!(geo.IsTransformed()))
			{
				for (int i = 0; i < 3; i++)
				{
					ofs	<< geo.CmNode().coords[i] / 1000 << ' ';
				}
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					ofs	<< matrix[i] << ' ';
				}
			}
			ofs	<< 0.0f << ' ' << 0.0f << std::endl 
				<< 0.0f << std::endl;
			for (int i = 0; i < 3; i++)
			{
				ofs	<< 0.0f << ' ';
			}
			ofs	<< 0.0f << ' ' << 0.0f << std::endl 
				<< 0.0f << std::endl;
			for (int i = 0; i < 9; i++)
			{
				ofs << matrix[i+3]
				<< (((i + 1) % 5 && (i != 8)) ? ' ' : '\n');
			}
		}
		ofs.close();
	}
	else
	{
		exceptions::ThrowFileNotOpened(FILE_DEFAULT_INITIAL_CONDITIONS);
	}
}

// static
void FModel::FullClean()
{
	Clean();

	//fs::DeleteFiles("*.gmr");
	fs::DeleteFiles("*.mov");
	fs::DeleteFiles("*.lba");
	fs::DeleteFiles("*.obj");
	fs::DeleteFiles("*.exe");
	fs::DeleteFiles("*.f31");
	fs::DeleteFiles("*.for");
	fs::DeleteFiles("*.F");
	fs::DeleteFiles("*.F32");
	fs::DeleteFiles("*.FOR");
	fs::DeleteFiles("*.lbm");
	fs::DeleteFiles("*.krt");
	fs::DeleteFiles("*.inc");
	fs::DeleteFiles("*.lst");
	fs::DeleteFiles("*.prg");
}

// static
void FModel::Clean()
{
	fs::DeleteFiles("*.sts");
	fs::DeleteFiles("*.fi");
	fs::DeleteFiles("*.f");
	fs::DeleteFiles("*.obj");
	fs::DeleteFiles("*.o");
	fs::DeleteFiles("*.F");
	fs::DeleteFiles("*.F32");
	fs::DeleteFiles("*.FOR");
	fs::DeleteFiles("*.mhl");
	fs::DeleteFiles("*.rhl");
	fs::DeleteFiles("*.mdl");
	fs::DeleteFiles("*.ras");
	fs::DeleteFiles(FILE_PMODEL_SPRINGM);
	fs::DeleteFiles(FILE_PMODEL_DEFAULMUPR);
	//TODO: lbm, krt, ...
}

bool IsRoad(FBody body)
{
	return body.Number() == 49 ? true : false;
}

/** Проверка налиия дороги в модели
* @return true если дорога есть
*/
bool FModel::IsRoadExist() const
{
	return find_if(_bodies.begin(), _bodies.end(), &IsRoad)!=_bodies.end();
}

/** Выбрать набор управляющих параметров для расчета
* @param paramsSetId - номер набора параметров
*/
void FModel::SetSolveControlParamsSet(int paramsSetId)
{
	_header.SetSolveParamsSet(paramsSetId);
}

/** Выбрать набор управляющих параметров для построения графиков
* @param paramsSetId - номер набора параметров
*/
void FModel::SetPostprocessorControlParamsSet(int paramsSetId)
{
	_header.SetPostprocessorControlParamsSet(paramsSetId);
}

void FModel::OutputAdditionalVehicleParams(ofstream& stream) const
{
	int numberOfLeftAxles = 0;
	int numberOfRightAxles = 0;
	int maxNumberOfLeftAxles=0;
	int maxNumberOfRightAxles=0;
	int isRoadHasLargeMoves = false;

	for(int i = 0; i < _bodies.size(); i++)
	{
		if(IsRoad(_bodies[i]))
		{
			FBodyDof bodyDof;
			_bodies[i].GetBodyDof(bodyDof);
			isRoadHasLargeMoves = !bodyDof.IsSmallMovementsMode();
		}
	}

	for(int i = 0; i < _joints.size(); i++)
	{
		int nodeNumber = -1;
		if(_joints[i].BodyNumber1() == 49)
			nodeNumber = _joints[i].NodeNumber1();
		else if(_joints[i].BodyNumber2() == 49)
			nodeNumber = _joints[i].NodeNumber2();
		if(nodeNumber != -1)
		{
			if(nodeNumber > 10)
			{
				numberOfRightAxles++;
				if(nodeNumber >	maxNumberOfRightAxles)
					maxNumberOfRightAxles = nodeNumber;
			}
			else
			{
				numberOfLeftAxles++;
				if(nodeNumber >	maxNumberOfLeftAxles)
					maxNumberOfLeftAxles = nodeNumber;
			}
		}

	}

	if(numberOfLeftAxles!=numberOfRightAxles)
		exceptions::ThrowMessage("Unequal number of left and right tracks contact points with road");
	//Разное количество точек контакта с дорогой левой и правой колей
	if(maxNumberOfLeftAxles > numberOfLeftAxles)
		exceptions::ThrowMessage("Error in ordering of left track connection points");
	//Нарушена последовательность перечисления узлов дороги левой колеи
	if(maxNumberOfRightAxles-10 > numberOfRightAxles)
		exceptions::ThrowMessage("Error in ordering of right track connection points");
	if(maxNumberOfLeftAxles && !isRoadHasLargeMoves)
		std::cout << "For small moves of road in BL36 distanses between axles not generated\n";
	//exceptions::ThrowMessage("For small moves of road in BL36 distanses between axles not generated");
	// При кинематическом возмущении для малых движений в блоке 36 не генерируются расстояния между мостами
	if(maxNumberOfLeftAxles)
	{
		stream << blockHeaders[10];
		float dlt=0;
		int l = 8;
		int number = 49;
		int id = 1;
		//"BL36 NPDS     NUZL    NPDS2    NUZL2     SKOL\n",			//10		
		for(int i=0; i<maxNumberOfLeftAxles; i++)
		{
			stream << 
				setw(9) << number << ' ' <<
				setw(l) << id	  << ' ' <<
				setw(l) << number << ' ' << 
				setw(l) << i+1	  << ' ' << 
				setw(l) << dlt	  << std::endl;
		}

		stream << blockHeaders[11];
		//"BL37 NKOL       HS        D       FZ       D1       D2\n"	//11
		//unlink("default.rhl");
		unlink("exlist.dat"); // TODO: выяснить что это
		//FHeader* header = this->Header();

		const double (&leftTrackParams)[6] = _header.GetSolveSpecialParams()->leftTrackParams;

		stream << 
			setw(9) << 1 << ' ' <<
			setw(l) << leftTrackParams[0] << ' ' <<
			setw(l) << leftTrackParams[1] << ' ' << 
			setw(l) << leftTrackParams[2] << ' ' << 
			setw(l) << leftTrackParams[3] << ' ' << 
			setw(l) << leftTrackParams[4] << std::endl;

		const double (&rightTrackParams)[6] = _header.GetSolveSpecialParams()->rightTrackParams;

		stream << 
			setw(l) << 2 << ' ' <<
			setw(l) << rightTrackParams[0] << ' ' <<
			setw(l) << rightTrackParams[1] << ' ' << 
			setw(l) << rightTrackParams[2] << ' ' << 
			setw(l) << rightTrackParams[3] << ' ' << 
			setw(l) << rightTrackParams[4] << std::endl;
	}
}

/**
* Загрузить начальные условия в геометрию тел
* @param icoParams - начальные условия
*/
/*
void FModel::LoadIco
	(
		IcoParams* icoParams
	)
{
	ifstream ifs(icoParams->icoFilePath);

	if (ifs.is_open())
	{
		int bodyNumber;
		int bodyIndex;
		float tmp;

		ifs	>> bodyNumber;
		while (!ifs.eof())
		{
			bodyIndex = GetIndexFromNumber(bodyNumber);
			if (bodyIndex != -1)
			{
				FGeometry* geo;

				geo = GetBodyGeometry(bodyIndex);

				float buf[12];

				for (int i = 0; i < 3; i++)
				{
					ifs	>> buf[i];
				}
				for (int i = 0; i < 9; i++)
				{
					ifs	>> tmp;
				}
				for (int i = 0; i < 9; i++)
				{
					ifs	>> buf[i + 3];
				}
				geo->SetTransformMtx(buf);
			}
			else
			{
				// TODO: выяснить, а нужно ли это здесь вообще?
				//exceptions::ThrowBodyNotFoundException(bodyNumber);

				for (int i = 0; i < 12; i++)
				{
					ifs	>> tmp;
				}
				for (int i = 0; i < 9; i++)
				{
					ifs	>> tmp;
				}
			}
			ifs	>> bodyNumber;
		}

		ifs.close();
	}
	else
	{
		exceptions::ThrowFileNotOpened(icoParams->icoFilePath);
	}
}
*/

/**
* Загрузить начальные условия в геометрию тел
* @param icoParams - начальные условия
*/
void FModel::LoadIco
	(
		const IcoParams* icoParams
	)
{
	MatrixMap::const_iterator it = icoParams->begin();

	while (it != icoParams->end())
	{
		int bodyId = GetBodyIndexByNumber(it->first);

		if (bodyId != -1)
		{
			FGeometry geo = GetBodyGeometry(bodyId);

			geo.SetTransformationMatrix(it->second);
			SetBodyGeometry(bodyId, geo);
		}
		++it;
	}
}

/** Выбрать набор управляющих параметров для расчета
* @param solveParamsSetId - номер набора параметров
*/
void FModel::SelectSolveParamsSet(int solveParamsSetId)
{
	if (solveParamsSetId != -1)
	{
		_header.SetSolveParamsSet(solveParamsSetId);
	}

	IcoParams* icoParams;

	icoParams = _header.GetIcoParams();
	LoadIco(icoParams);
}

void FModel::SelectVariablesGroup(int variablesSetId)
{
	if (variablesSetId != -1)
		_header.SetVariablesParamsSet(variablesSetId);
	Eval();
}

Calculator* FModel::GetCalculator() const
{
	Calculator* calc = new Calculator();
	FVariablesGroup fvar = _header.GetCurrentVariablesGroup();
	for (int i = 0; i < fvar.Size(); i++)
		calc->AddParameter(fvar[i].name.c_str(), fvar[i].expression.c_str());
	return calc;
}

SolverParamsBase* FModel::AddNewSolver(size_t bodyId, SolverTypes solverType)
{
	if(Header().GetMph() == nullptr)
		Header().AddMph();
	return GetBody(bodyId).AddNewSolver(solverType);
}

SolverParamsBase* FModel::GetSolver(size_t bodyId, size_t solverId)
{
	return GetBody(bodyId).GetMph()->GetSolver(solverId);
}

SolverIntParams& FModel::AddNewIntParams(SolverParamsBase* solver)
{
	solver->SetIntParamsId(GetMphHeader()->CountIntParamsSets());
	return GetMphHeader()->AddNewIntParams();
}

SolverGridParams& FModel::AddNewGridParams(SolverParamsBase* solver)
{
	solver->SetGridParamsId(GetMphHeader()->CountGridParamsSets());
	return GetMphHeader()->AddNewGridParams();
}

SolverSpecialParams& FModel::AddNewSpecialParams(SolverParamsBase* solver)
{
	solver->SetSpecialParamsId(GetMphHeader()->CountSpecialParamsSets());
	return GetMphHeader()->AddNewSpecialParams(solver->GetSolverType());
}

FMphHeader* FModel::AddMphHeader()
{
	if(Header().GetMph() != nullptr)
		delete Header().GetMph();
	Header().AddMph();
	return Header().GetMph();
}

/** Добавить новый набор параметров для решателя
* @param mapper - BcMapper
* @param bcType - тип граничных условий (от него зависит состав параметров)
* @param surfaceId - номер набора индексов поверхностей
*/
BcParams& FModel::AddNewBcParams(BcMapper& mapper, 
	BcTypes bcType, int surfaceId)
{
	if(surfaceId < 0)
	{
		mapper.SetRestBcId(GetMphHeader()->CountBcParamsSets());
		mapper.EnableBc();
	}
	else
	{
		mapper.SetBcParamId(surfaceId,GetMphHeader()->CountBcParamsSets());
	}
	
	BcParams& bcParams = GetMphHeader()->AddNewBcParams(bcType);
	bcParams.CreateDefault(bcType);

	
	return bcParams;
}

const FBody& FModel::GetBodyByNumber
	(
		int bodyNumber
	)	const
{
	auto bodyIterator = std::find_if
		(
			_bodies.begin(),
			_bodies.end(),
			[bodyNumber](const FBody& body)
			{
				return body.Number() == bodyNumber;
			}
		);

	if (bodyIterator == _bodies.end())
	{
		exceptions::ThrowNoElementInCollection("body", "bodies");
	}

	return *bodyIterator;
}

FBody& FModel::GetBodyByNumber
	(
		int bodyNumber
	)
{
	return const_cast<FBody&>(static_cast<const FModel*>(this)->GetBodyByNumber(bodyNumber));
}

bool FModel::IsBodyExistByNumber
	(
		int bodyNumber
	)	const
{
	auto bodyIterator = std::find_if
		(
			_bodies.begin(),
			_bodies.end(),
			[bodyNumber](const FBody& body)
			{
				return body.Number() == bodyNumber;
			}
		);

	return bodyIterator != _bodies.end();
}

const FForce& FModel::GetForceByNumber
	(
		int forceNumber
	)	const
{
	auto forceIterator = std::find_if
		(
			_forces.begin(),
			_forces.end(),
			[forceNumber](const FForce& force)
			{
				return force.Number() == forceNumber;
			}
		);

	if (forceIterator == _forces.end())
	{
		exceptions::ThrowNoElementInCollection("force", "forces");
	}

	return _forces[0];
}

FForce& FModel::GetForceByNumber
	(
		int forceNumber
	)
{
	return const_cast<FForce&>(static_cast<const FModel*>(this)->GetForceByNumber(forceNumber));
}

bool FModel::IsForceExistByNumber
	(
		int forceNumber
	)	const
{
	auto forceIterator = std::find_if
		(
			_forces.begin(),
			_forces.end(),
			[forceNumber](const FForce& force)
			{
				return force.Number() == forceNumber;
			}
		);

	return forceIterator != _forces.end();
}

const FJoint& FModel::GetJointByNumber
	(
		int jointNumber
	)	const
{
	auto jointIterator = std::find_if
		(
			_joints.begin(),
			_joints.end(),
			[jointNumber](const FJoint& joint)
			{
				return joint.Number() == jointNumber;
			}
		);

	if (jointIterator == _joints.end())
	{
		exceptions::ThrowNoElementInCollection("joint", "joints");
	}

	return _joints[0];
}

FJoint& FModel::GetJointByNumber
	(
		int jointNumber
	)
{
	return const_cast<FJoint&>(static_cast<const FModel*>(this)->GetJointByNumber(jointNumber));
}

bool FModel::IsJointExistByNumber
	(
		int jointNumber
	)	const
{
	auto jointIterator = std::find_if
		(
			_joints.begin(),
			_joints.end(),
			[jointNumber](const FJoint& joint)
			{
				return joint.Number() == jointNumber;
			}
		);

	return jointIterator != _joints.end();
}

FForceChar& FModel::GetForceCharByNumber( int forceCharNumber )
{
	for (int i=0; i< _forceChars.size(); i++)
	{
		if(_forceChars[i].Number() == forceCharNumber)
			return _forceChars[i];
	}
	exceptions::ThrowNoElementInCollection("forceChar","forceChars");
	return _forceChars[0];
}

const FForceChar& FModel::GetForceCharByNumber(int forceCharNumber) const
{
	for (int i=0; i< _forceChars.size(); i++)
	{
		if(_forceChars[i].Number() == forceCharNumber)
			return _forceChars[i];
	}
	exceptions::ThrowNoElementInCollection("forceChar","forceChars");
	return _forceChars[0];
}

const FJointChar& FModel::GetJointCharByNumber(int jointNumber) const
{
	for (int i=0; i < _jointChars.size(); i++)
	{
		if (_jointChars[i].Number() == jointNumber)
			return _jointChars[i];
	}
	exceptions::ThrowNoElementInCollection("jointChar", "jointChars");

	return _jointChars[0];
}

TransformationMatrix FModel::GetBodyTransformation(int bodyId) const
{
	return GetBodyGeometry(bodyId).GetTransformationMatrix();
}

const FGeometryPoint& FModel::GetBodyCm(int bodyId) const
{
	return GetBodyGeometry(bodyId).CmNode();
}

map<int, FMph*> FModel::GetBodyMphs()
{
	map<int, FMph*> mphs;

	for (size_t i = 0; i < _bodies.size(); i++)
	{
		FMph *pMph = _bodies[i].GetMph();

		if (pMph != nullptr)
			mphs[i] = pMph;
	}

	return mphs;
}

/**
* Добавить в карту "занятых" узлов очередной узел
* @param occupiedNodeNumbersMap - карта "занятых" узлов
* @param bodyNumber - номер тела, которому принадлежит узел
* @perem nodeNumber - номер добавляемого "занятого" узла
*/
void AddOccupiedNode
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap,
		int bodyNumber,
		int nodeNumber
	)
{
	BodyNodeNumbersMap::iterator bodyIterator;

	bodyIterator = occupiedNodeNumbersMap.find(bodyNumber);
	if (bodyIterator == occupiedNodeNumbersMap.end())
	{
		NodeNumbersSet nodeNumbers;

		nodeNumbers.insert(nodeNumber);
		occupiedNodeNumbersMap.insert
			(
				BodyNodeNumbersPair
					(
						bodyNumber,
						nodeNumbers
					)
			);
	}
	else
	{
		NodeNumbersSet& nodeNumbers = bodyIterator->second;
		
		if (nodeNumbers.find(nodeNumber) == nodeNumbers.end())
		{
			nodeNumbers.insert(nodeNumber);
		}
	}
}

/**
* Добавить в карту "занятых" узлов узлы, полученные от тел
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void FModel::AddOccupiedNodesFromBodies
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)
{
	for (int bodyIndex = 0; bodyIndex < _bodies.size(); bodyIndex++)
	{
		AddOccupiedNode
			(
				occupiedNodeNumbersMap,
				_bodies[bodyIndex].Number(),
				GetBodyGeometry(bodyIndex).CmNodeNumber() + 1
			);
	}
}

/**
* Добавить в карту "занятых" узлов узлы, полученные от сил
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void FModel::AddOccupiedNodesFromForces
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)	const
{
	for (int forceIndex = 0; forceIndex < _forces.size(); forceIndex++)
	{
		const FForce& force = _forces[forceIndex];

		AddOccupiedNode
			(
				occupiedNodeNumbersMap,
				force.BodyNumber(),
				force.NodeNumber()
			);
	}
}

/**
* Добавить в карту "занятых" узлов узлы, полученные от соединительных элементов
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void FModel::AddOccupiedNodesFromJoints
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)	const
{
	for (int jointIndex = 0; jointIndex < _joints.size(); jointIndex++)
	{
		const FJoint& joint = _joints[jointIndex];
		
		AddOccupiedNode
			(
				occupiedNodeNumbersMap,
				joint.BodyNumber1(),
				joint.NodeNumber1()
			);
		AddOccupiedNode
			(
				occupiedNodeNumbersMap,
				joint.BodyNumber2(),
				joint.NodeNumber2()
			);
	}
}

/* Добавить в карту "занятых" узлов узлы, полученные от макроса,
* в котором указаны номера узлов для некоторого тела
* @param macroSpecification - описание макроса
* @param modelParams - модельные параметры макроса
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void AddOccupiedNodesFromMacro
	(
		const FMacroSpecification& macroSpecification,
		const vector<string>& modelParams,
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)
{
	int paramIndex = 0;
	const MacroModelParam* modelParamSpecification;

	modelParamSpecification = macroSpecification.GetModelParam(0);
	do
	{
		while
			(
				(modelParamSpecification != nullptr) &&
				(!modelParamSpecification->IsBodyNumber())
			)
		{
			paramIndex++;
			modelParamSpecification = macroSpecification.GetModelParam(paramIndex);
		}
		if (modelParamSpecification == nullptr)
		{
			continue;
		}

		int bodyNumber;
		
		bodyNumber = atoi(modelParams.at(paramIndex).c_str());
		paramIndex++;
		modelParamSpecification = macroSpecification.GetModelParam(paramIndex);
		if
			(
				(modelParamSpecification == nullptr) ||
				(!modelParamSpecification->IsNodeNumber())
			)
		{
			continue;
		}

		do
		{
			int nodeNumber = atoi(modelParams.at(paramIndex).c_str());

			AddOccupiedNode
				(
					occupiedNodeNumbersMap,
					bodyNumber,
					nodeNumber
				);
			paramIndex++;
			modelParamSpecification = macroSpecification.GetModelParam(paramIndex);
		}
		while
			(
				(modelParamSpecification != nullptr) &&
				(modelParamSpecification->IsNodeNumber())
			);
	}
	while (modelParamSpecification != nullptr);
}

/* Добавить в карту "занятых" узлов узлы, полученные от специального макроса,
* который неявно использует номера узлов
* @param macroSpecification - описание макроса
* @param modelParams - модельные параметры макроса
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void AddOccupiedNodesFromSpecialMacro
	(
		const FMacroSpecification& macroSpecification,
		const vector<string>& modelParams,
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)
{
	if (macroSpecification.GetType() != 809) // tugol1.mac
	{
		return;
	}

	const MacroModelParam* modelParamSpecification;

	modelParamSpecification = macroSpecification.GetModelParam(0);
	for
		(
			int paramIndex = 1;
			modelParamSpecification != nullptr;
			modelParamSpecification = macroSpecification.GetModelParam(paramIndex++)
		)
	{
		if (!modelParamSpecification->IsBodyNumber())
		{
			continue;
		}
		if (modelParamSpecification->GetVariableName().find("Nbody") == string::npos)
		{
			continue;
		}

		int bodyNumber;
		
		bodyNumber = atoi(modelParams.at(paramIndex).c_str());
		AddOccupiedNode
			(
				occupiedNodeNumbersMap,
				bodyNumber,
				1
			);
	}
}

/**
* Добавить в карту "занятых" узлов узлы, полученные от макросов
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void FModel::AddOccupiedNodesFromMacros
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)
{
	MacroSpecificationMap macroSpecificationMap;

	if (!FromResourcesFile(macroSpecificationMap))
	{
		std::cout << "Errors may occur in the processing cut geometry!" << std::endl;

		return;
	}
	for (int macroIndex = 0; macroIndex < _header.GetMacroCount(); macroIndex++)
	{
		const FMacro& macro = _header.GetFMacro(macroIndex);
		MacroSpecificationMap::const_iterator macroIterator;
		
		macroIterator = macroSpecificationMap.find(macro.Type());
		if (macroIterator == macroSpecificationMap.end())
		{
			continue;
		}

		const FMacroSpecification& macroSpecification = macroIterator->second;
		const vector<string>& modelParams = macro.GetMParams();

		if (modelParams.size() != macroSpecification.GetModelParamsCount())
		{
			continue;
		}

		AddOccupiedNodesFromMacro
			(
				macroSpecification,
				modelParams,
				occupiedNodeNumbersMap
			);
		AddOccupiedNodesFromSpecialMacro
			(
				macroSpecification,
				modelParams,
				occupiedNodeNumbersMap
			);
	}
}

/**
* Найти номера узлов модели, к которым приложены силы или соединительные элементы
* @param occupiedNodeNumbersMap - карта "занятых" узлов
*/
void FModel::FindOccupiedNodes
	(
		BodyNodeNumbersMap& occupiedNodeNumbersMap
	)
{
	occupiedNodeNumbersMap.clear();
	AddOccupiedNodesFromBodies(occupiedNodeNumbersMap);
	AddOccupiedNodesFromForces(occupiedNodeNumbersMap);
	AddOccupiedNodesFromJoints(occupiedNodeNumbersMap);
	AddOccupiedNodesFromMacros(occupiedNodeNumbersMap);
}

void FModel::AddCurrentIcoParamsSet(const string& name)
{
	IcoParams icoParams;
	icoParams.name = name;
	icoParams.icoFilePath = name + EXT_ICO;
	for(int i=0; i < _bodies.size(); i++)
	{
		float mtx[12];
		GetBodyGeometry(i).GetTransformArray(mtx);
		TransformationMatrix matrix(mtx);
		icoParams.InsertIcoBodyMatrix(_bodies[i].Number(), matrix);
	}
	icoParams.Save();
	Header().AddCurrentIcoParamsSet(name, icoParams);
}

FGeometryPoint FModel::GetLocalTransformedNode(size_t bodyId, size_t nodeId) const
{
	return GetBodyGeometry(bodyId).GetLocalTransformedNode(nodeId);
}

const FGeometryPoint& FModel::GetTransformedNode(const BodyNodeNumber& bodyNode) const
{
	int bodyId = GetBodyByNumber(bodyNode.bodyNumber).Id();
	return GetBodyGeometry(bodyId).GetNode(bodyNode.nodeNumber-1); // id is zero-based
}

int FModel::AddBody(const FBody& body)
{
	_bodies.push_back(body);
	_numbersToIndexes[body.Number()] = body.Id();
	return body.Id();
}

int FModel::AddJoint(const FJoint& joint)
{
	int id = _joints.size();
	_joints.push_back(joint);
	return id;
}

int FModel::AddJointChar(const FJointChar& jointChar)
{
	int id = _jointChars.size();
	_jointChars.push_back(jointChar);
	return id;
}

int FModel::AddForce(const FForce& force)
{
	int id = _forces.size();
	_forces.push_back(force);
	return id;
}

int FModel::AddForceChar(const FForceChar& forceChar)
{
	int id = _forceChars.size();
	_forceChars.push_back(forceChar);
	return id;
}

int FModel::AddGeometry(const FGeometry& geometry)
{
	int id = _geometries.size();
	_geometries.push_back(geometry);
	return id;
}

void FModel::AddMph(const FMph& mph)
{
	_mphs.push_back(mph);
}

void FModel::SetBody(const FBody& body, int id)
{
	_bodies[id] = body;
}

void FModel::SetGeometry(const FGeometry& geometry, int id)
{
	_geometries[id] = geometry;
}

/**
* Найти имя файла с моделью в текущей директории
* (подразумевается, что в одной директории м.б. только одна модель и ее backup)
* @param pathToModel - путь к директории с моделью (если пусто, то используется текущая рабочая директория)
* @return найденное имя файла с моделью (если модель не найдена, возвращается пустая строка)
*/
// static
string FModel::FindModelFileName
	(
		const string& pathToModel
	)
{
	string path;

	if (pathToModel.empty())
	{
		fs::DirectoryLister dl(".", "*.frm");
		string curPath;

		dl.Init();
		while (!dl.IsLast())
		{
			dl.Next(curPath);
			if (curPath.find("backup") == string::npos)
			{
				path = curPath;
				break;
			}
		}
	}
	else
	{
		path = pathToModel;
	}

	return path;
}

/**
* Загрузить модель, расположенную в данной директории
* @param pathToModel - путь к директории с моделью (если пусто, то используется текущая рабочая директория)
* @return указатель на загруженную модель (если модель не загружена, выбрасывается исключение)
*/
// static
FModel* FModel::LoadModel
	(
		const string& pathToModel
	)
{
	string modelName = FindModelFileName(pathToModel);

	return new FModel(modelName.c_str());
}


void FModel::Eval() const
{
	Calculator* calc = GetCalculator();

	for (int i = 0; i < _bodies.size(); i++)
		_bodies[i].Eval(*calc);
	for (int i = 0; i < _joints.size(); i++)
		_joints[i].Eval(*calc);
	for (int i = 0; i < _jointChars.size(); i++)
		_jointChars[i].Eval(*calc);

	delete calc;
}

