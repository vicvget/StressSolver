#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <vector>
#include <cstdio>
#include <string>

using std::vector;
using std::string;


/*
Класс содержащий граничные условия для сетки описанной 
объектом OccRectilinearGrid
*/
class BoundaryCondition
{
private:
	string myName;
	string myVar;
	string myId;
	vector<int> myPoints;
public:
	int operator [] (int id) const {return myPoints[id];}
	size_t CountPoints() const {return myPoints.size();}

	BoundaryCondition(void);
	BoundaryCondition(string name, string var, string id, const vector<int>& points):
		myName(name),
		myVar(var),
		myId(id),
		myPoints(points)
		{};
	~BoundaryCondition(void);
	/*
	* @return возвращет значение граничного условия 
	*/
	const string& GetVar() const {return myVar;}
	void SetVar(string s){myVar=s;}
	/*
	* @return возвращет идентификатор граничного условия 
	*/
	const string& GetID() const {return myId;}
	void SetID(string s){myId=s;}
	/*
	* @return возвращет имя граничного условия 
	*/
	const string& GetName() const {return myName;}
	void SetName(std::string nm) {myName = nm;}
	const vector<int>& GetPoints() const {return myPoints;}
	void SetPoints(vector<int> x){ myPoints = x;}
	void AddPoint(int i){myPoints.push_back(i);}
	int GetPoint(size_t i){return myPoints.at(i);}
	size_t GetNumBP() const {return myPoints.size();}
};

// список граничных условий для сетки, описанной объектом класса OccRectilinearGrid
typedef vector<BoundaryCondition*> BoundaryConditions;

#endif