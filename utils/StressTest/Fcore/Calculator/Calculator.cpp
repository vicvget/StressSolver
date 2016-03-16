#include "Calculator.h"

#include "../../fcore/wrappers/fcUnistd.h"
#include "../../fcore/wrappers/StringRoutines.h"

#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>


namespace Calc
{

Calculator::Calculator()
	:
		input_string(nullptr),
		cur_pos(),
		current_token(UNKNOWN),
		number_value()
{
	std::fill(name_string, name_string + sizeof(name_string), 0);
	param p;
	strcpy(p.name,"pi");
	p.expression[0] = 0;
	p.value=M_PI;
	params_table.push_back(p);
	strcpy(p.name,"e");
	p.expression[0] = 0;
	p.value=2.7182818284590452354;
	params_table.push_back(p);
	strcpy(p.name,"g");
	p.expression[0] = 0;
	p.value=9.81;
	params_table.push_back(p);
	sp.kpt = 0;
	sp.kpu = 500;
	sp.kmx = 0;
	sp.kse = 0;
	sp.kss = 0;
	sp.kspring = 0;
}

// Добавить системный параметр
double Calculator::AddSystemParam(const char *s)
{
	double value = 0;
	if(!strnicmp(s,"npt",3))
	{
		sp.kpt++;
		if(sp.kpt==49) sp.kpt++;
		value = sp.kpt;
	}
	else if(!strnicmp(s,"npu",3)){
		value = ++sp.kpu;
	}
	else if(!strnicmp(s,"nmx",3)){
		value = ++sp.kmx;
	}
	else if(!strnicmp(s,"nse",3)){
		value = ++sp.kse;
	}
	else if(!strnicmp(s,"nss",3)){
		value = ++sp.kss;
	}
	else
	{
		throw CEParamNotFound(s, cur_pos);
	}
	AddParameter(s,value);
	return value;
}

/** Вычислить выражение
*
* @param expression - вычисляемое выражение
* @return вычисленное значение
*/
double Calculator::Eval(const char* expression)
{
	if (strlen(expression) == 0)
	{
		throw CEParamNotFound("empty");
	}

	input_string = expression;
	cur_pos = 0;
	get_token();
	if(current_token == END)
		return 0;
	else
		return expr();
}

double Calculator::TryEval(const char* expression)
{
	double result = 0.;
	try
	{
		result = Eval(expression);
	}
	catch (CalculatorException* ex)
	{
		std::cerr << "EVAL EXCEPTION:" << ex->ToString() << std::endl;
	}
	return result;
}

/** Вычислить выражение
*
* @param expression - вычисляемое выражение
* @param value - вычисленное значение
* @return true, если вычислено успешно
*/
bool Calculator::TryEval(const char* expression, double& value)
{
	try
	{
		value = Eval(expression);
		return true;
	}
	catch (const CalculatorException&)
	{
		return false;
	}
}

bool Calculator::Check(const char* expression)
{
	double value;
	return TryEval(expression, value);
}

// Добавить параметр-константу
void Calculator::AddParameter(const char* name, double value)
{
	param* par = look(name);
	if(par == nullptr)
	{
		par = new param;
		strcpy(par->name, name);
		par->expression[0] = '\0';
		par->value = value;
		params_table.push_back(*par);
		delete par;

	}
	else
	{
		par->value = value;
		UpdateAfter(name);		
	}

}

// Добавить параметр-выражение
void Calculator::AddParameter(const char* name, const char* expression)
{
	param* par = look(name);
	if(par == nullptr)
	{
		par = new param;
		strcpy(par->name, name);
		strcpy(par->expression, expression);
		par->value = Eval(expression);
		params_table.push_back(*par);
		delete par;

	}
	else
	{
		strcpy(par->expression, expression);
		UpdateAfter(name);		
	}
}

// Пересчитать параметры
void Calculator::Update()
{
	vector<param>::iterator li = this->params_table.begin();
	while (li != params_table.end())
	{
		if (li->expression[0] != '\0')
			li->value = Eval(li->expression);
		++li;
	}
}

// Пересчитать параметры до заданного включительно
void Calculator::UpdateTo(const char* name)
{
	vector<param>::iterator li = this->params_table.begin();
	while(li != params_table.end())
	{
		if(li->expression[0] != '\0')
			li->value = Eval(li->expression);
		if(!strcmp(li->name, name))
			break;
		++li;
	}
}

// Пересчитать параметры после заданного включительно
void Calculator::UpdateAfter(const char* name)
{
	vector<param>::iterator li = this->params_table.begin();
	int eval = 0;
	while(li != params_table.end())
	{
		if((li->expression[0] != '\0')&&(eval))
			li->value = Eval(li->expression);
		if(!strcmp(li->name, name))
			eval = 1;
		++li;
	}
}
// Поиск параметра в таблице по имени
param* Calculator::look(const char* name)
{
	vector<param>::iterator li = params_table.begin();
	while(li != params_table.end())
	{
		if(!strcmp(li->name, name))
			return &(*li);
		++li;
	}
	return nullptr;
}

// Устанавливает следующий TOKEN и возвращает его
int Calculator::get_token()
{
	char ch;
	do
	{
		ch = input_string[cur_pos++];
		if(ch == '\0')
			return END;		
	}
	while(ch != '\n' && (ch == isspace(ch)|| ch == ' '));

	switch(ch)
	{
	case '!':
		return current_token = END;
	case '\n':
	case ';':
		cur_pos++;
		return current_token = PRINT;
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	case '0':
	case '.':
		{
			char* p = name_string;
			while(isdigit(ch) || ch == 'e' || ch == 'E' || ch == '.')
			{
				*p++ = ch;
				if((ch == 'e')||(ch == 'E'))
				{
					ch = input_string[cur_pos++];
					if((ch == '-')||(ch == '+'))
					{
						*p++ = ch;
						ch = input_string[cur_pos++];
					}
				}
				else 
					ch = input_string[cur_pos++];
			}
			cur_pos--;
			*p = 0;

			int k = sscanf(name_string, "%lf", &number_value);

			if(k == 0)
			{
				throw CENotANumber(name_string, cur_pos);
			}

			return current_token = NUMBER;
		}
	case '+':
	case '-':
	case '/':
	case '*':
	case '%':
	case '^':
	case '(':
	case ')':
	case '=':
		return current_token = token_value(ch);
	default:
		// 2010 FIX for XYZ bug with variables _px
		if(isalpha(ch)||ch=='_')
		{
			char* p = name_string;
			//			ch =input_string[cur_pos++];
			while(isalnum(ch)||ch=='_')
			{
				*p++ = ch;
				ch = input_string[cur_pos++];
			}
			cur_pos--;
			*p = 0;
			return current_token = NAME;
		}
		else
		{
			std::cout <<"name str " <<name_string << std::endl;
			throw CENotAName(name_string, cur_pos);
		}
	}	
}

// Вычислитель уровня 0
double Calculator::prim()
{
	switch(current_token)
	{
	case NUMBER:
		get_token();
		return number_value;
		break;
	case NAME:
		{
			struct param* par;

			if (get_token() == ASSIGN)
			{
				int isRequireAdd = 0;

				par = look(name_string);
				if(par == nullptr)
				{
					par = new param;
					strcpy(par->name,name_string);
					isRequireAdd = 1;
				}
				strcpy(par->expression,input_string+cur_pos);
				get_token();
				double result = expr();
				par->value=result;
				if(isRequireAdd)
				{
					params_table.push_back(*par);
					delete par;
				}
				return result;
			}
			else
			{
			if(current_token == LP)
			{
				if(!strcmp(name_string,"cos"))
					return cos(prim());
				else if(!strcmp(name_string,"sin"))
					return sin(prim());
				else if(!strcmp(name_string,"exp"))
					return exp(prim());
				else if(!strcmp(name_string,"tan"))
				{
					double tmp = prim();
					if((tmp - floor(tmp/M_PI)*M_PI)<FLOATZERO)
					{
						throw CEInvalidValue("Incorrect param of TAN", tmp, cur_pos);
					}
					return tan(tmp);
				}
				else if(!strcmp(name_string,"sqrt"))
				{
					double tmp = prim();
					if(tmp < 0)
					{
						throw CEInvalidValue("Negative param of SQRT", tmp, cur_pos);
					}
					return sqrt(tmp);
				}
				else if(!strcmp(name_string,"ln"))
				{
					double tmp = prim();
					if(tmp < FLOATZERO)
					{
						throw CEInvalidValue("Non-positive param of LN", tmp, cur_pos);
					}
					return log(tmp);
				}
				else
				{
					throw CEInvalidFunction(name_string, cur_pos);
				}
			}
			else
			{
				par = look(name_string);
				if(par != nullptr)
				{
					return par->value;
				}
				else
				{
					return AddSystemParam(name_string);
					//throw CEParamNotFound(name_string, cur_pos);
					// AddSysParr
				}
				return 0;
			}
			}
		}
	case MIN:
		{
			get_token();
			return -prim();
		}
		break;
	case LP:
		{
			get_token();
			double e=expr();
			if(current_token != RP)
			{
				throw CENoBracket(cur_pos);
			}
			get_token();
			return e;
		}
	case END:
		return 1;
	default:
		throw CEInvalidToken(current_token, cur_pos);
	}
}

// Вычислитель уровня 1
double Calculator::tterm()
{
	double left = prim();
	while(1)
	{
		switch (current_token)
		{
		case POW:
			get_token();			
			left=pow(left, prim());			
			break;
		default:
			return left;
		}
	}
}


// Вычислитель уровня 2
double Calculator::term()
{
	double left = tterm();
	while(1)
	{
		switch (current_token)
		{
		case MUL:
			get_token();
			left*=tterm();
			break;
		case DIV:
			{
				get_token();
				double tmp = tterm();
				if(tmp == 0)
				{
					throw CEDivisionByZero(cur_pos);
				}
				left /= tmp;
				break;
			}
		case RES:
			{
				get_token();
				double tmp = tterm();
				if(tmp == 0)
				{
					throw CEDivisionByZero(cur_pos);
				}
				left = left - floor(left/tmp)*tmp;
				break;
			}
		default:
			return left;
		}
	}
}


// Вычислитель уровня 3
double Calculator::expr()
{
	double left = term();
	while(1)
	{
		switch (current_token)
		{
		case ADD:
			get_token();
			left+=term();
			break;
		case MIN:
			get_token();
			left-=term();
			break;
		default:
			return left;
		}
	}
}


// Отладочный листинг таблицы параметров
void Calculator::Print()
{
	vector<param>::iterator li = params_table.begin();
	while(li != params_table.end())
	{
		std::cout << li->name << '\t' << li->value << '\t' << li->expression << '\n';
		++li;
	}
}

//Проверить наличие параметра
double Calculator::GetParameter(const char* name)
{
	param* p = look(name);

	if (p == nullptr)
	{
		throw CEParamNotFound(name);
	}

	return p->value;
}

//Проверить наличие параметра
int Calculator::CheckParameter(const char* name, double* value)
{
	param* p = look(name);
	if(p == nullptr)
		return 0;
	if(value != nullptr)
		*value = p->value;
	return 1;
}

char* Calculator::GetLastParameter()
{
	if(params_table.size() == 0)
	{
		return nullptr;
	}
	else
	{
		vector<param>::iterator li = params_table.end();
		--li;
		return li->name;
	}
}


// Удалить добавленные после заданного параметры
void Calculator::RemoveAfter(const char* name)
{
	vector<param>::iterator li = params_table.begin();
	while(li != params_table.end())
	{
		if(!strcmp(li->name, name))
		{
			++li;
			break;
		}
		++li;
	}
	if(li != params_table.end())
	{
		params_table.erase(li, params_table.end());
	}
}

// Удалить параметр
void Calculator::RemoveParam(const char* name)
{
	vector<param>::iterator li = params_table.begin();
	while(li != params_table.end())
	{
		if(!strcmp(li->name, name))
		{
			break;
		}
		++li;
	}
	if(li != params_table.end())
	{
		params_table.erase(li);
	}
}

void Calculator::Push()
{
	vector<param> state;
	state.assign(params_table.begin(),params_table.end());
	//params_table.cop
	params_states.push(state);
}

void Calculator::Pop()
{
	if(!params_states.empty())
	{
		vector<param> state = params_states.top();
		params_table.assign(state.begin(), state.end());
		params_states.pop();
	}
}

}