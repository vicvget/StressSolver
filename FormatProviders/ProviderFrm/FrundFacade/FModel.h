#ifndef FMODEL_H

#define FMODEL_H


#include "FHeader.h"
#include "FMphHeader.h"
#include "FBody.h"
#include "FJoint.h"
#include "FJointChar.h"
#include "FForce.h"
#include "FCunit.h"
#include "FForceChar.h"
#include "FFreeSolver.h"
#include "MphParamsSet\FreeSolverGridParams.h"

using Calc::Calculator;


// карта конверсии одних индексов в другие
// ключ - конвертируемый индекс
// значение - сконвертированный индекс
typedef pair<int, int> IndexConversionPair;
typedef map<int, int> IndexConversionMap;

// список номеров узлов
typedef set<int> NodeNumbersSet;

// карта номеров узлов, соответствующих определенным телам
// ключ - номер тела
// значение - список номеров узлов, принадлежащих данному телу
typedef pair<int, NodeNumbersSet> BodyNodeNumbersPair;
typedef map<int, NodeNumbersSet> BodyNodeNumbersMap;


/**
* Класс для считывания модели ФРУНД
*
* @author Getmanskiy Victor
*/
class FModel
{
private:
	// заголовок
	FHeader _header;
	//заголовок междисциплинарных парарметров
	//FMphHeader fMphHeader;
	// механические характеристики связей
	vector<FJointChar> _jointChars;
	// связи между телами
	vector<FJoint> _joints;
	// механические характеристики сил
	vector<FForceChar> _forceChars;
	// внешние силы, действующие на тела системы
	vector<FForce> _forces;
	// тела
	vector<FBody> _bodies;
	// управляющие элементы (связь с макросами)
	vector<FCunit> _cunits;
	// геометрия
	vector<FGeometry> _geometries;
	// мультифизика
	vector<FMph> _mphs;

	// название (для совместимости не должно превышать 8 символов?)
	string _name;
	//путь к каталогу
	string _path;
	// карта конверсии номеров тел в их индексы в массиве bodies
	IndexConversionMap _numbersToIndexes;

	//свободные решатели
	vector<FFreeSolver> _freeSolvers;

	FDebugParams _debugParams;

	/**
	* Формирование карты конверсии номеров тел в индексы в массиве bodies
	*/
	void FormNumbersToIndexes();

	/**
	* Добавить в карту "занятых" узлов узлы, полученные от тел
	* @param occupiedNodeNumbersMap - карта "занятых" узлов
	*/
	void AddOccupiedNodesFromBodies
		(
			BodyNodeNumbersMap& occupiedNodeNumbersMap
		);

	/**
	* Добавить в карту "занятых" узлов узлы, полученные от сил
	* @param occupiedNodeNumbersMap - карта "занятых" узлов
	*/
	void AddOccupiedNodesFromForces
		(
			BodyNodeNumbersMap& occupiedNodeNumbersMap
		)	const;

	/**
	* Добавить в карту "занятых" узлов узлы, полученные от соединительных элементов
	* @param occupiedNodeNumbersMap - карта "занятых" узлов
	*/
	void AddOccupiedNodesFromJoints
		(
			BodyNodeNumbersMap& occupiedNodeNumbersMap
		)	const;

	/**
	* Добавить в карту "занятых" узлов узлы, полученные от макросов
	* @param occupiedNodeNumbersMap - карта "занятых" узлов
	*/
	void AddOccupiedNodesFromMacros
		(
			BodyNodeNumbersMap& occupiedNodeNumbersMap
		);

public:
	FModel();

	void Load(const char* path);

	FModel(const char* path);

	/** Cохранение по заданному пути
	* Создает файлы модели в формате ФРУНДа: frm, и .elf по fileId
	* @param path - путь к файлу
	*/
	void Save(const string& path) const;

	/** Cохранение в текущей директории под имененм в _name
	* Создает файлы модели в формате ФРУНДа: frm, и .elf по fileId
	*/
	void Save() const;

	/** Проверка наличия дороги в модели
	* @return true если дорога есть
	*/
	bool IsRoadExist() const;


#pragma region getters

	/**
	* Получить наименование файла модели с расширением
	* @return наименование файла модели с расширением
	*/
	string GetModelFileName() const;

	/** Возвращает количество тел в модели
	* @return количество тел в модели
	*/
	size_t CountBodies() const { return _bodies.size(); }

	/** Возвращает количество тел в модели
	* @return количество тел в модели
	*/
	size_t GetBodiesCount() const { return _bodies.size(); }

	/** Возвращает количество геометрий в модели
	* @return количество геометрий в модели
	*/
	size_t GetGeometriesCount() const { return _geometries.size(); }

	/** Возвращает количество сил в модели
	* @return количество тел в модели
	*/
	size_t GetForcesCount() const { return _forces.size(); }

	/** Возвращает количество характеристик сил в модели
	* @return количество характеристик сил в модели
	*/
	size_t GetForceCharsCount() const { return _forceChars.size(); }

	/** Возвращает количество сединительных элементов
	* @return количество соединительных элементов
	*/
	size_t GetJointsCount() const {return _joints.size();}

	/** Возвращает количество характеристик сединительных элементов
	* @return количество характеристик соединительных элементов
	*/
	size_t GetJointCharsCount() const {return _jointChars.size();}

	/** Возвращает количество cunit
	* @return количество cunit
	*/
	size_t GetCunitsCount() const {return _cunits.size();}

	/** Возвращает количество свободных решателей
	* @return количество свободных решателей
	*/
	size_t GetFreeSolversCount() const { return _freeSolvers.size(); }

	/**
	* Возвращает индекс тела в массиве bodies по его номеру
	* @param bodyNumber - номер тела
	* @return индекс тела в массиве bodies (-1, если тело с таким номером не найдено)
	*/
	int GetBodyIndexByNumber
		(
			int bodyNumber
		)	const;

	/** Получить параметры тела
	* @param bodyId - индекс тела
	* @return параметры тела
	*/
	const FBody& GetBody(int bodyId) const {return _bodies[bodyId];}

	FBody& GetBody(int bodyId) {return _bodies[bodyId];}

	const FBody& GetBodyByNumber
		(
			int bodyNumber
		)	const;

	FBody& GetBodyByNumber
		(
			int bodyNumber
		);

	bool IsBodyExistByNumber
		(
			int bodyNumber
		)	const;

	FGeometry& GetGeometry(int geometryId) {return _geometries[geometryId];}
	const FGeometry& GetGeometry(int geometryId) const {return _geometries[geometryId];}

	/** Получить параметры cunit
	* @param cunitId - индекс cunit
	* @return параметры cunit
	*/
	FCunit& GetCunit(int cunitId) {return _cunits[cunitId];}

	/** Получить параметры соединительного элемента
	* @param jointId - индекс соединительного элемента
	* @return параметры соединительного элемента
	*/
	FJoint& GetJoint(int jointId) {return _joints[jointId];}
	const FJoint& GetJoint(int jointId) const {return _joints[jointId];}
	
	const FJoint& GetJointByNumber
		(
			int jointNumber
		)	const;
	
	FJoint& GetJointByNumber
		(
			int jointNumber
		);

	bool IsJointExistByNumber
		(
			int jointNumber
		)	const;

	/** Получить параметры силы
	* @param forceId - индекс силы
	* @return параметры силы
	*/
	FForce& GetForce(int forceId) {return _forces[forceId];}
	const FForce& GetForce(int forceId) const {return _forces[forceId];}

	const FForce& GetForceByNumber
		(
			int forceNumber
		)	const;

	FForce& GetForceByNumber
		(
			int forceNumber
		);

	bool IsForceExistByNumber
		(
			int forceNumber
		)	const;

	/** Получить параметры характеристики соединительного элемента
	* @param jointId - индекс характеристики соединительного элемента
	* @return параметры характеристики соединительного элемента
	*/
	const FJointChar& GetJointChar(int jointCharId) const {return _jointChars[jointCharId];}
	FJointChar& GetJointChar(int jointCharId) { return _jointChars[jointCharId]; }

	const FJointChar& GetJointCharByJoint(const FJoint& joint) const
	{
		return GetJointCharByNumber(joint.CharNumber());
	}

	const FJointChar& GetJointCharByNumber(int jointNumber) const;
	
	/** Получить параметры характеристики силы
	* @param forceId - индекс характеристики силы
	* @return параметры характеристики силы
	*/
	FForceChar& GetForceChar(int forceCharId) {return _forceChars[forceCharId];}
	const FForceChar& GetForceChar(int forceCharId) const { return _forceChars[forceCharId]; }
	
	const FForceChar& GetForceChar(const FForce& force) const { return GetForceCharByNumber(force.CharNumber()); }
	const FForceChar& GetForceCharByNumber(int forceCharNumber) const;
	FForceChar& GetForceCharByNumber(int forceCharNumber);

	/** Получить геометрию для тела (sc)
	* @param bodyId - индекс тела
	* @return геометрия
	*/
	const FGeometry& GetBodyGeometry(int bodyId) const 
	{ 
		return _geometries.at(_bodies.at(bodyId).GeometryId());
	}

	TransformationMatrix GetBodyTransformation(int bodyId) const;

	const FGeometryPoint& GetBodyCm(int bodyId) const;

	void SetBodyGeometry(int bodyId, const FGeometry& geometry)
	{
		_geometries[_bodies[bodyId].GeometryId()] = geometry;
	}

	/** Получить параметры для междисциплинарного расчета для тела
	* @param bodyId - индекс тела
	* @return параметры для междисциплинарного расчета
	*/
	//FMultiphysics* GetBodyMultiphysics(int bodyId) { return bodies[bodyId].GetMultiphysics();}


	/** Получить параметры для междисциплинарного расчета для тела
	* @param bodyId - индекс тела
	* @return параметры для междисциплинарного расчета
	*/
	FMph* GetBodyMph(int bodyId) { return _bodies[bodyId].GetMph();}

	/** Получить параметры для междисциплинарного расчета для тел
	* @return параметры для междисциплинарного расчета
	*/
	map<int, FMph*> GetBodyMphs();

	// TODO:
	//FMph* GetBodyMph(int bodyId) { return bodies[bodyId].GetMph();}

	/** Получить параметры для междисциплинарного расчета для свободного решателя
	* @return параметры для междисциплинарного расчета для свободного решателя
	*/
	FFreeSolver& GetFreeSolver(int freeSolverId) { return _freeSolvers[freeSolverId]; }

#pragma endregion

	/** Возвращвет вектор из 6 компонент, true, если есть данная степень свободы
	* первые 3 поступательные, последние три - вращательные.
	* поступательные изображаются как стрелки, вращательные - как моменты (дуга вокруг оси со стрелкой)
	*
	* Например, если есть сила по x и z и момент вокруг Oy, то будет true,false,true,false,true,false
	* @param components - компоненты силы
	* @param forceId - Id силы
	*/
	void GetComponentsVector(vector<bool>& components, int forceId) const
	{
		int charNumber = _forces[forceId].CharNumber();

		if ((charNumber > 0) && (charNumber <= static_cast<int>(_forceChars.size())))
		{
			_forceChars[charNumber - 1].GetComponentsVector(components);
		}
		else
		{
			// TODO: разобраться с этой проблемой !!! Почему CharNumber() может быть > chars.Size()!
			//exceptions::ThrowCollectionOutOfBounds("ComponentVector");
		}
	}
	void GetForceComponentsVector(vector<bool>& components,int forceNumber) 
	{
		int charNumber = GetForceByNumber(forceNumber).CharNumber();

		GetForceCharByNumber(charNumber).GetComponentsVector(components);
	}

	/** Вывод расчетной схемы в формате ФРУНД
	* @param stream - выходной поток
	*/
	void Output(ofstream& stream, int stage);

	/** Вывод в текстовый файловый поток заголовка расчетной схемы
	* @param stream - выходной поток
	*/
	void OutputHeader(ofstream& stream, int stage) const;
	
	/** Вывод в текстовый файловый поток переменных текущего набора
	* @param stream - выходной поток
	*/
	void OutputVariables(ofstream &stream, int stage) const;
	
	/** Вывод в текстовый файловый поток макросов
	* @param stream - выходной поток
	*/
	void OutputMacros(ofstream &stream, int stage) const;
	
	/** Вывод информации о телах в модели
	* @param stream - выходной поток
	*/
	void OutputBodies(ofstream& stream, int stage) const;
	
	/** Вывод информации о фиктивных телах в модели
	* @param stream - выходной поток
	*/
	void OutputCunits(ofstream& stream, int stage) const;
	
	/** Вывод информации о характеристиках шарниров в модели
	* @param stream - выходной поток
	*/
	void OutputJointChars(ofstream& stream, int stage) const;
	
	/** Вывод информации о шарнирах в модели
	* @param stream - выходной поток
	*/
	void OutputJoints(ofstream& stream, int stage) const;
	
	/** Вывод информации о характеристиках сил в модели
	* @param stream - выходной поток
	*/
	void OutputForceChars(ofstream& stream, int stage) const;
	
	/** Вывод информации о силах в модели
	* @param stream - выходной поток
	*/
	void OutputForces(ofstream& stream, int stage) const;
	
	/** Вывод информации о силах в модели (BL36-37)
	* @param stream - выходной поток
	*/
	void OutputAdditionalVehicleParams(ofstream& stream) const;

	/** Вывод управляющих параметров для расчета(UPRF)
	* @param stream - выходной поток
	*/
	void OutputSolveControlParams(ofstream &stream) const;

	/** Вывод дополнительных управляющих параметров
	* @param stream - выходной поток
	*/
	void OutputSolveExtraControlParams(ofstream &stream) const;

	/** Вывод параметров моделей с дорогами
	*/
	void OutputSolveParams() const;
	
	/** Вывод начальных условий, соответствующих исходной геометрии
	*/
	void OutputDefaultIco() const;
	
	/** Выбрать набор управляющих параметров для расчета
	* @param paramsSetId - номер набора параметров
	*/
	void SetSolveControlParamsSet(int paramsSetId);

	/** Выбрать набор управляющих параметров для построения графиков
	* @param paramsSetId - номер набора параметров
	*/
	void SetPostprocessorControlParamsSet(int paramsSetId);

	/**
	* Загрузить начальные условия в геометрию тел
	* @param icoParams - начальные условия
	*/
	void LoadIco
		(
			const IcoParams* icoParams
		);

	/** Выбрать набор управляющих параметров для расчета
	* @param solveParamsSetId - номер набора параметров
	*/
	void SelectSolveParamsSet(int solveParamsSetId = -1);

	/** Обновить численые данные в модели в соответствии с набором переменных
	* @param variablesSetId - номер набора переменных, если -1, то использовать текущий
	*/	
	void SelectVariablesGroup(int variablesSetId = -1);
	
	/** Вывод управляющих параметров для анализа (UPRF)
	* @param ofs - выходной поток
	*/
	void OutputAnazControlParams(ofstream &ofs) const;
	
	/** Вывод информации о телах в формате библиотеки ФРУНД (lba)
	*/
	void SaveLba() const;
	
	// Вывод информации о телах в формате библиотеки ФРУНД (mdl)
	//	void SaveLbm(AvmodelInterface& avmodel, char* path);

	// Сформировать последовательность вывода в sequence
	// type = 0 - для тел, sequence = bodySequence[rest]
	// type = 1 - для тел, по кол-ву уравнений связей (bl[0]=-1)
	// type = 2 - для тел, по кол-ву уравнений связей (bl[0]=-2)
	// type = 3 - для с.е., в обратном порядке следования тел bodySequence
	void SortOutputSequence(int type);

	/**
	* Загрузить геометрию профиля из *949.gmr или *49.gmr файлов
	* @return признак успешной (true) или неуспешной (false) загрузки
	*/
	bool LoadRoadProfileGeometry();

#pragma region partition
	/** Вывод информации о модели в виде графа в формате совместимом с GraphViz
	* http://www.graphviz.org/
	*/
	void OutputModelToGraph(ofstream& stream);


	/** Возвращает модель в формате CSR 
	* http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf
	*/
	void GetGraphCSR(int *xadj ,int *adjncy,int *vwgt);

	/** Индексирование модели
	*/
	void IndexModel();
#pragma endregion

#pragma region appenders
	int AddBody(const FBody& body);
	int AddJoint(const FJoint& joint);
	int AddJointChar(const FJointChar& jointChar);
	int AddForce(const FForce& force);
	int AddForceChar(const FForceChar& forceChar);
	int AddGeometry(const FGeometry& geometry);
	void AddMph(const FMph& mph);
#pragma endregion

#pragma region modifiers
	void SetBody(const FBody& body, int id);
	void SetGeometry(const FGeometry& geometry, int id);
#pragma endregion

	
#pragma region accessors
	string Name() const { return _name; }
	void Name(string val) { _name = val; }

	string Path() const { return _path; }
	void Path(string val) { _path = val; }

	vector<FBody> Bodies() const { return _bodies; }
	void Bodies(vector<FBody> val) { _bodies = val; }

	vector<FJoint> Joints() const { return _joints; }
	void Joints(vector<FJoint> val) { _joints = val; }

	vector<FForce> Forces() const { return _forces; }
	void Forces(vector<FForce> val) { _forces = val; }

	vector<FForceChar> ForceChars() const { return _forceChars; }
	void ForceChars(vector<FForceChar> val) { _forceChars = val; }

	vector<FJointChar> JointChars() const { return _jointChars; }
	void JointChars(vector<FJointChar> val) { _jointChars = val; }
	
	vector<FCunit> Cunits() const { return _cunits; }

	vector<FFreeSolver> FreeSolvers() const { return _freeSolvers; }
	void FreeSolvers(vector<FFreeSolver> val) { _freeSolvers = val; }

#pragma endregion

	/** Очистка каталога модели
	*/
	static
	void Clean();

	/** Очистка каталога модели в том числе и от старых файлов
	*/
	static
	void FullClean();

	/** Создание калькулятора с текущим набором переменных
	* @return калькулятор
	*/
	Calculator* GetCalculator() const;

	/** Получить заголовок
	* @return заголовок
	*/
	const FHeader& Header() const {return _header;}

	/** Получить заголовок
	* @return заголовок
	*/
	FHeader& Header() {return _header;}

	/** Получить междисциплинарный заголовок
	*/
	FMphHeader *GetMphHeader() {return _header.GetMph();}

	/*Добавить текущий набор паметров как начальные условия
	*/
	void AddCurrentIcoParamsSet(const string& name);

	/** Добавить рештель для тела
	* @param bodyId - номер тела
	* @param solverType - тип решателя
	* @return решатель
	*/
	SolverParamsBase* AddNewSolver(size_t bodyId, SolverTypes solverType);
	SolverParamsBase* GetSolver(size_t bodyId, size_t solverId);
	
	SolverIntParams& AddNewIntParams(SolverParamsBase* solver);
	SolverGridParams& AddNewGridParams(SolverParamsBase* solver);
	SolverSpecialParams& AddNewSpecialParams(SolverParamsBase* solver);
	
	/** Добавить новый набор параметров для решателя
	* @param mapper - BcMapper
	* @param bcType - тип граничных условий (от него зависит состав параметров)
	* @param surfaceId - номер набора индексов поверхностей, -1 = restBcParams
	*/
	BcParams& AddNewBcParams(BcMapper& mapper, 
		BcTypes bcType, int surfaceId=-1);

	/** Добавить и вернуть новый заголовок междисциплинарных параметров
	* @return заголовок междисциплинарных параметров
	*/
	FMphHeader* AddMphHeader();

	/** Возвращает узел с учетом локальной трансформации (для начальных условий)
	* @param bodyId - номер тела
	* @ nodeId - номер узла в transformedNodes
	*/
	FGeometryPoint GetLocalTransformedNode(size_t bodyId, size_t nodeId) const;

	/** Возвращает узел с учетом локальной трансформации (для начальных условий)
	* @param bodyId - номер тела
	* @ nodeId - номер узла в transformedNodes
	*/
	const FGeometryPoint& GetTransformedNode(const BodyNodeNumber& bodyNode) const;
	
	/** Сохранение в указанный каталог
	* @param newPath - каталог
	*/
	void SaveTo(const string& newPath) const;

	/**
	* Найти номера узлов модели, к которым приложены силы или соединительные элементы
	* @param occupiedNodeNumbersMap - карта "занятых" узлов
	*/
	void FindOccupiedNodes
		(
			BodyNodeNumbersMap& occupiedNodeNumbersMap
		);
	

	/**
	* Найти имя файла с моделью в данной директории
	* (подразумевается, что в одной директории м.б. только одна модель и ее backup)
	* @param pathToModel - путь к директории с моделью (если пусто, то используется текущая рабочая директория)
	* @return найденное имя файла с моделью (если модель не найдена, возвращается пустая строка)
	*/
	static
	string FindModelFileName
		(
			const string& pathToModel = {}
		);

	/**
	* Загрузить модель, расположенную в данной директории
	* @param pathToModel - путь к директории с моделью (если пусто, то используется текущая рабочая директория)
	* @return указатель на загруженную модель (если модель не загружена, выбрасывается исключение)
	*/
	static
	FModel* LoadModel
		(
			const string& pathToModel = {}
		);

	/**
	* Вычисляет значения констант в элементах модели
	*/
	void Eval() const;

};


#endif // FMODEL_H
