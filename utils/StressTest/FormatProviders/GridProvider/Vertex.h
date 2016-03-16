#ifndef GRAPHPART_H

#define GRAPHPART_H

//#pragma GCC visibility push(hidden)
#include <fstream>
#include <iostream>
#include <vector>


using std::vector;


struct VertexPoint
{

	double x{};

	double y{};

	double z{};

};

class Vertex
{
public:

	vector<int> adjVert;
	vector<int> adjWeight;
	int _localNumber{};
	int _globalNumber{};
	int size{};
	double weight{};
	VertexPoint pt;

	Vertex()
	{
		adjVert.reserve(6);
	}

	Vertex
		(
			const Vertex* other
		)
		:
			// TODO: разобраться, почему не копируются adjWeight и weight
			adjVert(other->adjVert),
			_localNumber(other->_localNumber),
			_globalNumber(other->_globalNumber),
			size(other->size),
			pt(other->pt)
	{
	}

	bool operator == (const Vertex &other) const
	{
		bool flag = true;

		if (_globalNumber == other._globalNumber && size == other.size && pt.x == other.pt.x && pt.y == other.pt.y && pt.z == other.pt.z)
		{
			for (size_t i = 0; i < adjVert.size(); i++)
			{
				if (adjVert[i] != other.adjVert[i])
				{
					flag = false;
				}
			}
		}
		else
		{
			flag = false;
		}
		return flag;
	}

	bool operator != (const Vertex &other) const
	{
		return !(*this == other);
	}

};

//#pragma GCC visibility pop
#endif