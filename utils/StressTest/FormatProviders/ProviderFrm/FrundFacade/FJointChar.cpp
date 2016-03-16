#include "FJointChar.h"

#include "../../../fcore/fcore.h"
#include "../../../AdditionalModules/Workarounds/GNU_gcc/string/to_string.h"

#include <functional>


using std::setw;


FJointChar::FJointChar()
	:
		_maxStiffness(MaxStiffness)
{
	
}

FJointChar::FJointChar
	(
		int orderIndex
	)
	:
		FElement(orderIndex),
		_maxStiffness(MaxStiffness)
{
}

bool FJointChar::Load(ifstream& stream)
{
	fs::ReadLineString(stream, _name);
	stream >> _number;

	int directions[6], nDirections = 0, isDirection;

	for (int i = 0; i < 6; i++)
	{
		stream >> isDirection;
		if(isDirection)
		{
			directions[nDirections++] = i + 1;
		}
	}
	for (int i = 0; i < nDirections; i++)
	{
		for (int j = 0; j < 2; j++)
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
			j == 0 ? _springParams.push_back(component) : _dampingParams.push_back(component);
		}
	}
	_maxStiffness = MaxStiffness;

	return true;
}

void FJointChar::Eval(Calc::Calculator &calc) const
{
	for(int i = 0; i < _springParams.size(); i++)
	{
		for(int j = 0; j < 2; j++)
		{
			const CharComponent& comp = j == 0 ? _springParams[i] : _dampingParams[i];
			if(j == 0)
			{
				comp.stiffness = calc.TryEval(comp.params[0].c_str());
			}
		}
	}
}

void FJointChar::Save(ofstream& ofs) const
{
	ofs << _name << std::endl
		<< _number << std::endl;
	
	int directions[6];
	memset(directions, 0, sizeof(int)*6);
	int i, j;
	for(i = 0; i < _springParams.size(); i++)
	{
		directions[_springParams[i].direction-1] = _springParams[i].type;		
	}
	for(i = 0; i < 6; i++)
	{
		ofs << directions[i] << std::endl;
	}
	for(i = 0; i < _springParams.size(); i++)
	{
		ofs << _springParams[i].type << ' '
			<< _springParams[i].params.size() << std::endl;
		for(j = 0; j < _springParams[i].params.size(); j++)
			ofs << _springParams[i].params[j] << std::endl;	
		ofs << _dampingParams[i].type << ' '
			<< _dampingParams[i].params.size() << std::endl;
		for(j = 0; j < _dampingParams[i].params.size(); j++)
			ofs << _dampingParams[i].params[j] << std::endl;	
	}
}

void FJointChar::Output
	(
		ofstream& stream,
		int isRType
	)	const
{
	if (isRType)
	{
		// BL33 NMX    PNAPR    PSLAG     KPAR
		const int l = 8;

		stream << "@ \'" << _name << "\'\n";
		for (int i = 0; i < _springParams.size(); i++)
		{
			for (int j = 0; j < 2; j++)
			{
				const CharComponent& comp = j == 0 ? _springParams[i] : _dampingParams[i];
				const vector<string>& compParams = comp.params;
				string stiffness = compParams[0];
				vector<std::reference_wrapper<const string>> componentParams(compParams.begin(), compParams.end());

				componentParams[0] = stiffness;
				if (j == 0)
				{
					// TODO: ???
					if (_maxStiffness < comp.stiffness)
					{
						stiffness = to_string(_maxStiffness);
					}
				}

				int addStr = (componentParams.size() % 4) ? 1 : 0;

				stream
					<< setw(11) << _number << ' '
					<< setw(l) << comp.direction << ' '
					// << setw(l) << '1'  << ' '
					<< setw(l) << j + 1 << ' '
					<< setw(l) << componentParams.size() / 4 + addStr << std::endl;
					// записывается не количество параметров, а количество строк по 4 под параметры!
					// << setw(l) << componentParams.size() << '\n';
				for (int k = 0; k < componentParams.size(); k++)
				{
					if (k % 4 == 0)
					{
						stream << setw(9) << componentParams[k].get();
					}
					else
					{
						stream << setw(l) << componentParams[k].get();
					}
					if (k % 4 == 3)
					{
						stream << '\n';
					}
					else
					{
						stream << ' ';
					}
				}

				int n = (4 - componentParams.size() % 4) % 4;

				for (int k = 0; k < n ; k++)
				{
					if (k != n - 1)
					{
						stream << setw(l) << "0" << ' ';
					}
					else
					{
						stream << setw(l) << "0" << '\n';
					}
				}
			}
		}
	}
	else
	{
		// BL16 NMX NAPR  TUF  KPAR  TDF  KPAR
		const int l = 4;

		stream << "@ \'" << _name << "\'\n";
		for(int i = 0; i < _springParams.size(); i++)
		{
			stream << setw(9) << _number << ' ' <<
				setw(l) << _springParams[i].direction << ' ' <<
				setw(l) << _springParams[i].type << ' ' << 
				setw(l) << _springParams[i].params.size() << ' ' <<
				setw(l) << _dampingParams[i].type << ' ' << 
				setw(l) << _dampingParams[i].params.size() << '\n';
		}
	}
}

void FJointChar::AddCharComponent(const CharComponent& charSpringComponent, const CharComponent& charDampingComponent)
{
	if(charSpringComponent.direction == charDampingComponent.direction)
	{
		_springParams.push_back(charSpringComponent);
		_dampingParams.push_back(charDampingComponent);
	}
	else
	{
		throw std::exception();
	}
}
