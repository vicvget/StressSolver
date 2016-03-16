#ifndef SURFACEPARAM_H

#define SURFACEPARAM_H


//#pragma GCC visibility push(hidden)
#include "Point3D.h"

#include <cstdio>
#include <vector>
#include <string>


using std::string;
using std::vector;


/**
* ����� ��� �������� ������ � �����������(�����, ����� ����������� � ���������� ��������� ����������)
*/
class SurfaceParam
{
protected:

	string name;
	string var;
	string id;

	//������ ����� ����������� �� �����������
	vector<Point3D *> SurfaceChildren;

	int index{};

public:

	SurfaceParam
		(
			const string& n,
			const string& v
		);

	void ClearPoints();
	//�������� ��������� �����
	void AddBoundaryPoints(int pntx,int pnty,int pntz);
	//������� ��������� ����� ��� �������� npnt
	int GetBoundaryPoints(int npnt,int *x,int *y,int *z);
	//������� ���������� ��������� �����
	size_t GetNumBP() const;

	/*
	* @return ��������� �������� ���������� ������� 
	*/
	const string& GetVar() const
	{
		return var;
	}

	void SetVar
		(
			const string& s
		)
	{
		var = s;
	}

	/*
	* @return ��������� ������������� ���������� ������� 
	*/
	const string& GetID() const
	{
		return id;
	}

	void SetID
		(
			const string& s
		)
	{
		id=s;
	}

	/*
	* @return ��������� ��� ���������� ������� 
	*/
	const string& GetName() const
	{
		return name;
	}

	void SetName
		(
			const string& nm
		)
	{
		name = nm;
	}

	const vector<Point3D *>& GetPoints() const
	{
		return SurfaceChildren;
	}

	void SetPoints
		(
			const vector<Point3D *>& x
		)
	{
		SurfaceChildren = x;
	}

	void SetIndex
		(
			int ind
		);

	int GetIndex() const;

};
//#pragma GCC visibility pop
#endif