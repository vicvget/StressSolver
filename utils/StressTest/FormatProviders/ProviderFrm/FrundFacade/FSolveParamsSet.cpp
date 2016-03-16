#include "FSolveParamsSet.h"

#include "../../../fcore/Exceptions/fcExceptions.h"
#include "../../../fcore/wrappers/FileRoutines.h"
#include "../../../fcore/fcore.h"

#include <sstream>
#include <iomanip>


using std::setw;
using std::endl;


#define XPLUS "Xplus"
#define XMINUS "Xminus"


void FSolveParamsSet::FullClean()
{
	icoParamsSets.clear();
	controlParamsSets.clear();
	extraControlParamsSets.clear();
	specialParamsSets.clear();
	paramSetIdsSets.clear();
}


FSolveParamsSet::FSolveParamsSet()
{
	icoParamsSets.resize(1);
	icoParamsSets[0].Init();
	controlParamsSets.resize(2);
	controlParamsSets[0].Init();
	controlParamsSets[1].Init(true);
	extraControlParamsSets.resize(2);
	extraControlParamsSets[0].Init();
	extraControlParamsSets[1].Init(true);
	specialParamsSets.resize(1);
	specialParamsSets[0].Init();
	paramSetIdsSets.resize(2);
	paramSetIdsSets[0].Init();
	paramSetIdsSets[1].Init(true);
}

FSolveParamsSet::~FSolveParamsSet()
{
}


void FSolveParamsSet::Read(ifstream& ifs)
{
	size_t count;
	
	fs::ReadLineValue(ifs, count);
	icoParamsSets.resize(count+1);
	vector<IcoParams>::iterator it = icoParamsSets.begin();
	it->Init();
	++it;
	if(it != icoParamsSets.end())
	{
		fs::ReadLineString(ifs, it->name);
		fs::ReadLineString(ifs, it->icoFilePath);
		it->Load();
		++it;
	}
	
	while(it != icoParamsSets.end())
	{
		fs::ReadLineString(ifs, it->name);
		fs::ReadLineString(ifs, it->icoFilePath);
		//Sergev & Timur =>
		it->Load();
		++it;
	}

	fs::ReadLineValue(ifs, count);
	controlParamsSets.resize(count+1);
	vector<ControlParams>::iterator itc = controlParamsSets.begin();
	itc->Init();
	++itc;

	int paramsCount = 0;
	while(itc != controlParamsSets.end())
	{
		fs::ReadLineString(ifs, itc->name);
		fs::ReadLineValue(ifs, paramsCount); // dummy
		fs::ReadLineValue(ifs, itc->numberOfSolveVariants);
		fs::ReadLineValue(ifs, itc->velocity);
		fs::ReadLineValue(ifs, itc->velocityStep);
		fs::ReadLineValue(ifs, itc->solverStep);
		fs::ReadLineValue(ifs, itc->solverOutputStep);
		fs::ReadLineValue(ifs, itc->solveType);
		// TODO: =8 what is it???
		//fs::ReadLineValue(ifs, itc->time); 
		fs::ReadLineValue(ifs, itc->time);
		fs::ReadLineValue(ifs, itc->outputFromTime);
		fs::ReadLineValue(ifs, itc->outputCode);
		++itc;
	}

	fs::ReadLineValue(ifs, count);
	extraControlParamsSets.resize(count+1);

	vector<ExtraControlParams>::iterator ite = extraControlParamsSets.begin();

	ite->Init();
	++ite;
	while (ite != extraControlParamsSets.end())
	{
		fs::ReadLineString(ifs, ite->name);
		//fs::ReadLineString(ifs, paramsCount); // dummy
		fs::ReadLineValue(ifs, paramsCount);
		fs::ReadLineValue(ifs, ite->factorOfNondissipativeRegulateForces);
		fs::ReadLineValue(ifs, ite->factorOfDissipativeRegulateForces);
		fs::ReadLineValue(ifs, ite->maxResidualInJoints);
		fs::ReadLineValue(ifs, ite->factorOfAMTFrequency);
		fs::ReadLineValue(ifs, ite->factorOfResidualRegulation);
		fs::ReadLineValue(ifs, ite->constraintsApplyStep);
		for (int i = 0; i < paramsCount - 6; i++)
		{
			int dummy;
			fs::ReadLineValue(ifs, dummy);
		}
		++ite;
	}

	fs::ReadLineValue(ifs, count);
	specialParamsSets.resize(count+1);

	vector<SpecialParams>::iterator its = specialParamsSets.begin();

	its->Init();
	++its;
	while (its != specialParamsSets.end())
	{
		fs::ReadLineString(ifs, its->name);
		fs::ReadLineValue(ifs, paramsCount); // dummy
		fs::ReadLineString(ifs, its->leftTrackFilePath);
		fs::ReadLineString(ifs, its->rightTrackFilePath);
		fs::ReadLineValue(ifs, its->useInitialVelocity);
		string strDirection;
		fs::ReadLineString(ifs, strDirection);
		its->direction = strDirection.compare(XPLUS) ? false : true;
		fs::ReadLineValue(ifs, its->trackType);
		for(int i = 0; i < 6; i++)
			fs::ReadLineValue(ifs, its->leftTrackParams[i]);
		for(int i = 0; i < 6; i++)
			fs::ReadLineValue(ifs, its->rightTrackParams[i]);
		fs::ReadLineValue(ifs, its->assumeAllJointsFlexible);
		++its;
	}

	fs::ReadLineValue(ifs, count);
	paramSetIdsSets.resize(count+1);

	vector<ParamSetIds>::iterator itp = paramSetIdsSets.begin();

	itp->Init();
	++itp;
	while (itp != paramSetIdsSets.end())
	{
		fs::ReadLineString(ifs, itp->name);
		stringstream line;
		fs::ReadLine(ifs, line);
		line >> itp->IcoParamsId >> 
			itp->ControlParamsId >>
			itp->ExtraControlParamsId >> 
			itp->SpecialParamsId >> 
			itp->VariableParamsId;
		itp->IcoParamsId--;
		itp->ControlParamsId--;
		itp->ExtraControlParamsId--;
		itp->SpecialParamsId--;
		itp->VariableParamsId--;
		++itp;
	}
	fs::ReadLineValue(ifs, currentId);
	if (currentId >= paramSetIdsSets.size())
	{
		currentId = paramSetIdsSets.size() - 1;
	}
}

void FSolveParamsSet::Write(ofstream& ofs) const
{
	// первый набор параметров "default" генерируется при считывании, его не надо записывать
	ofs << icoParamsSets.size()-1 << endl;
	vector<IcoParams>::const_iterator it = icoParamsSets.begin();
	// сохранем default.ico
	it->Save();
	++it;
	while (it != icoParamsSets.end())
	{
		ofs << it->name << endl << it->icoFilePath << endl;
		it->Save();
		++it;
	}
	ofs << controlParamsSets.size() - 1 << endl;

	vector<ControlParams>::const_iterator itc = controlParamsSets.begin();

	++itc;
	while (itc != controlParamsSets.end())
	{
		ofs << itc->name << endl <<
			'9' << endl << // params count
			itc->numberOfSolveVariants << endl <<
			itc->velocity << endl <<
			itc->velocityStep << endl <<
			itc->solverStep << endl <<
			itc->solverOutputStep << endl <<
			itc->solveType << endl <<
			itc->time << endl <<
			itc->outputFromTime << endl <<
			itc->outputCode << endl;
		++itc;
	}
	ofs << extraControlParamsSets.size()-1 << endl;

	vector<ExtraControlParams>::const_iterator ite = extraControlParamsSets.begin();

	++ite;
	while (ite != extraControlParamsSets.end())
	{
		ofs << ite->name << endl <<
			'6' << endl << // params count
			ite->factorOfNondissipativeRegulateForces << endl <<
			ite->factorOfDissipativeRegulateForces << endl <<
			ite->maxResidualInJoints << endl <<
			ite->factorOfAMTFrequency << endl <<
			ite->factorOfResidualRegulation << endl <<
			ite->constraintsApplyStep << endl;
		++ite;
	}
	ofs << specialParamsSets.size() - 1 << endl;

	vector<SpecialParams>::const_iterator its = specialParamsSets.begin();

	++its;
	while (its != specialParamsSets.end())
	{
		string strDirection = its->direction ? XPLUS : XMINUS;
		ofs << its->name << endl <<
			"18" << endl << // params count
			its->leftTrackFilePath << endl <<
			its->rightTrackFilePath << endl <<
			its->useInitialVelocity << endl <<
			strDirection << endl <<
			its->trackType << endl;
		for(int i = 0; i < 6; i++)
			ofs << its->leftTrackParams[i] << endl;
		for(int i = 0; i < 6; i++)
			ofs << its->rightTrackParams[i] << endl;
		ofs << its->assumeAllJointsFlexible << endl;
		++its;
	}
	ofs << paramSetIdsSets.size() - 1 << endl;

	vector<ParamSetIds>::const_iterator itp = paramSetIdsSets.begin();

	++itp;
	while (itp != paramSetIdsSets.end())
	{
		ofs << itp->name << endl <<
			itp->IcoParamsId+1 << ' ' <<
			itp->ControlParamsId+1 << ' ' <<
			itp->ExtraControlParamsId+1 << ' '<<
			itp->SpecialParamsId+1 << ' ' <<
			itp->VariableParamsId+1 << endl;
		++itp;
	}
	ofs << currentId << endl;
}


/** Вывод в uprf Параметров расчета
* @param ofs - output file stream
*/
void FSolveParamsSet::OutputControlParams(ofstream& ofs) const
{
	const ControlParams& controlParams = controlParamsSets[paramSetIdsSets[currentId].ControlParamsId];

	std::cout << MSG_OUT_CONTROL << std::endl
	<< "Control Params Group Name:" << controlParams.name << std::endl
	<< "Solve Type:" << controlParams.solveType << std::endl
	<< "Id of case:" << currentId << std::endl
	<< "Id of control:" << paramSetIdsSets[currentId].ControlParamsId << std::endl;
	// TODO: с этим разобраться
	
	ofs << setw(4) << controlParams.numberOfSolveVariants << ' ';
	ofs << setw(2) << controlParams.velocity << ' ';
	ofs << setw(4) << controlParams.velocityStep << ' ';
	ofs << setw(6) << controlParams.solverStep << ' ';
	ofs << setw(7) << controlParams.solverOutputStep << ' ';
	ofs << setw(6) << controlParams.solveType << ' ';
	ofs << setw(3) << controlParams.time*controlParams.velocity << ' ';
	ofs << setw(6) << controlParams.outputFromTime*controlParams.velocity << ' ';
	ofs << setw(7) << 1 << ' '; // TODO: ???
	ofs << setw(5) << controlParams.outputCode << ' ';
	ofs << "\nKVAR  V   DV      H      HZ   KODR   S     SO  KVARNU  KODN\n";
}


/** Вывод в epsilon.dat дополнительных отладочных параметров
* @param ofs - output file stream
*/
void FSolveParamsSet::OutputExtraParams(ofstream& ofs) const
{
	const ExtraControlParams& extraControlParams = 
		extraControlParamsSets[paramSetIdsSets[currentId].ExtraControlParamsId];
	ofs << extraControlParams.factorOfNondissipativeRegulateForces << endl
	 << extraControlParams.factorOfDissipativeRegulateForces << endl
	 << 0.0	<< endl	//?
	 << 1e-10 << endl //?
	 << extraControlParams.maxResidualInJoints << endl
	 << 1 << endl //?
	 << extraControlParams.factorOfResidualRegulation << endl
	 << extraControlParams.factorOfAMTFrequency << endl
	 << 2 << endl //?
	 << 250 << endl //?
	 << 0.5 << endl //?
 	 << 1 << endl //?
	 << extraControlParams.constraintsApplyStep << endl
	 << 1 << endl; //?

	//TODO: дописать комментарии на жестко 
	// забитые параметры ( в FSHELL так же)вывод параметров
}

/** Вывод в wat.cnt имен файлов с левым и правым профилем дороги
* @param ofs - outputt file stream
*/
void FSolveParamsSet::OutputRoadParams(ofstream& ofs) const
{
	const SpecialParams& specialParams = 
		specialParamsSets[paramSetIdsSets[currentId].SpecialParamsId];
	string leftTrackFilePath = fs::SplitFileFromPath(specialParams.leftTrackFilePath);
	string rightTrackFilePath = fs::SplitFileFromPath(specialParams.rightTrackFilePath);
	if(leftTrackFilePath != "none")
	{
		leftTrackFilePath = fs::CombinePathEnv(string(DIR_DATA), leftTrackFilePath);
		if(!fs::FileExist(leftTrackFilePath))
			exceptions::ThrowFileNotFound(leftTrackFilePath);
	}
	if(rightTrackFilePath != "none")
	{
		rightTrackFilePath = fs::CombinePathEnv(string(DIR_DATA), rightTrackFilePath);
		if(!fs::FileExist(rightTrackFilePath))
			exceptions::ThrowFileNotFound(rightTrackFilePath);
	}

	ofs << leftTrackFilePath << std::endl
	    << rightTrackFilePath << std::endl;
	// TODO: что делать с остальными параметрами ???
}

/** Вывод в default.ico начальных условий
* вывод в default.ico начальных условий
* @param ofs - output file stream
* @param id - номер ico
*/
void FSolveParamsSet::OutputIcoParams(ofstream& ofs, size_t id) const
{
	const IcoParams& icoParams = icoParamsSets[id];
	if((icoParams.icoFilePath != string(FILE_INITIAL_CONDITIONS)))
	{
		if(icoParams.icoFilePath == string(FILE_DEFAULT_INITIAL_CONDITIONS))
		{
			ofstream tmp(FILE_INITDEFA);
			tmp.close();
		}
		else
		{
			unlink(FILE_INITDEFA);
		}

		ifstream ifs(icoParams.icoFilePath);
		if(ifs.is_open())
		{
			string line;
			while(fs::ReadLineStringNoException(ifs, line))
			{
				ofs << line << endl;
			}
			ifs.close();
		}
		else
		{
			exceptions::ThrowFileNotOpened(icoParams.icoFilePath);
		}
	}
}

void FSolveParamsSet::OutputDefaultIcoParams() const
{
	ofstream ofs;
	ofs.open(FILE_INITIAL_CONDITIONS);
	if(ofs.is_open())
	{
		OutputIcoParams(ofs,0);
		ofs.close();
	}
	else
		exceptions::ThrowFileNotOpened(FILE_INITIAL_CONDITIONS);

}


/** Вывод в default.ico начальных условий
* вывод в default.ico начальных условий
* @param ofs - output file stream
*/
void FSolveParamsSet::OutputIcoParams(ofstream& ofs) const
{
	OutputIcoParams(ofs, paramSetIdsSets[currentId].IcoParamsId);
}