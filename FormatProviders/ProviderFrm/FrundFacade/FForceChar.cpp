#include "FForceChar.h"

#include "../../../fcore/fcore.h"


using std::setw;


FForceChar::FForceChar():FElement()
{
}

FForceChar::FForceChar(int orderIndex):FElement(orderIndex)
{
}

bool FForceChar::Load(ifstream& stream)
{
	fs::ReadLineString(stream, _name);
	stream >> _number;

	int directions[6], nDirections = 0, isDirection;

	for (int i = 0; i < 6; i++)
	{
		stream >> isDirection;
		if (isDirection)
		{
			directions[nDirections++] = i + 1;
		}
	}
	for (int i = 0; i < nDirections; i++)
	{
		CharComponent component;

		component.direction = directions[i];
		stream >> component.type;

		int nParams;
		string param;

		stream >> nParams;
		while (nParams > 0)
		{
			stream >> param;
			component.params.push_back(param);
			nParams--;
		}
		_params.push_back(component);
	}

	return true;
}

void FForceChar::Save(ofstream& ofs) const
{
	ofs << _name << std::endl;
	ofs << _number << std::endl;

	int directions[6];

	memset(directions, 0, 6 * sizeof(int));
	for (int i = 0 ; i < _params.size(); i++)
	{
		directions[_params[i].direction - 1] = _params[i].type;
	}
	for (int i = 0; i < 6; i++)
	{
		ofs << directions[i] << std::endl; //TODO: проверить 3
	}
	for (int i = 0 ; i < _params.size(); i++)
	{
		ofs << _params[i].type << ' '
			<< _params[i].params.size() << std::endl;
		for (int j = 0; j < _params[i].params.size(); j++)
		{
			ofs << _params[i].params[j] << std::endl;
		}
	}
}

void FForceChar::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	if (isRType) 
	{
		// BL34 NSS     NAPR     KSTR
		int l = 8;
		stream << "@ \'" << _name << "\'\n";
		for(int i = 0; i < _params.size(); i++)
		{
			const CharComponent& comp = _params[i];
			int addStr = (comp.params.size() % 3) ? 1 : 0;

			stream << setw(9) << _number << ' ' <<
				setw(l) << comp.direction << ' ' <<
				setw(l) << comp.params.size()/3+addStr << '\n';
			for(int k = 0; k < comp.params.size(); k++)
			{
				if( k % 3 == 0 )
					stream << setw(9) << comp.params[k];
				else
					stream << setw(l) << comp.params[k];
				if( k % 3 == 2 ) 
					stream << '\n';
				else 
					stream << ' ';
			}
			int n = (3 - comp.params.size() % 3) % 3;
			for(int k = 0; k < n ; k++)
			{
				if(k != n - 1) 
					stream << setw(l) << '0' << ' ';
				else  
					stream << setw(l) << '0' << '\n';

			}
		}
	}
	else
	{
		// BL18 NSS   NAPR    NTF   KPAR
		int l = 4;
		stream << "@ \'" << _name << "\'\n";
		for(int i = 0; i < _params.size(); i++)
		{
			stream << setw(9) << _number << ' ' <<
				setw(l) << _params[i].direction << ' ' <<
				setw(l) << _params[i].type << ' ' << 
				setw(l) << _params[i].params.size() << '\n';
		}
	}
}

void FForceChar::GetComponentsVector(vector<bool>& components) const
{
	components.clear();
	components.resize(6);
	for(int i = 0; i < _params.size(); i++)
	{
		components[i] = false;
	}
	for(int i = 0; i < _params.size(); i++)
	{
		components[_params[i].direction-1] = true;
	}
}
