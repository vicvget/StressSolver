#ifndef CALCULATOR_H

#define CALCULATOR_H


#ifdef GNUCPP
#	include "../../fcore/wrappers/fcUnistd.h"
#else
#	include <vector>
#	include <stack>
#	include <fstream>
#endif

#include "CalculatorExceptions.h"
#include "../../fcore/wrappers/fcUnistd.h"
#include "../../fcore/wrappers/FileRoutines.h"

#define FLOATZERO 1E-10

#ifndef M_PI
#	define M_PI 3.141592653589793238462643
#endif

//TODO: нормализовать документирование кода
//TODO: добавить arcsin, arccos, arctg, ...
namespace Calc
{
	/** Структура для хранения параметров
	*/
	struct param
	{
		// имя параметра
		char name[200];
		// выражение для параметра
		char expression[256];
		// значение параметра
		double value;	
	};

	/** Системные параметры калькулятора
	*/
	struct sysparams
	{
		int kpt; //
		int kpu; //
		int kmx; //количество механических характеристик;
		int kse; //
		int kss; //
		int kspring; //

		void write(ofstream& ofs) const
		{
			ofs << kpt << ' ' <<
				kpu << ' ' <<
				kmx << ' ' <<
				kse << ' ' <<
				kss << ' ' <<
				kspring << std::endl;
		}

		string getString() const
		{
			stringstream ss;
			ss << kpt << ' ' <<
				kpu << ' ' <<
				kmx << ' ' <<
				kse << ' ' <<
				kss << ' ' <<
				kspring;
			char buf[80];
			ss.getline(buf,80);
			string res = buf;
			return res;
		}

		void read(stringstream& ss)
		{
			ss >> kpt >> kpu >> kmx >> kse >> kss >> kspring;		
		}

		void read(const string& str)
		{
			stringstream ss;
			ss << str;
			read(ss);
		}

		void read(ifstream& ifs)
		{
			stringstream ss;
			fs::ReadLine(ifs, ss);
			read(ss);
		}
	};

	/** типы операций
	*/
	enum token_value 
	{
		UNKNOWN,
		NAME,
		NUMBER,
		PRINT,
		END,
		ADD = '+',
		MIN = '-',
		MUL = '*',
		DIV = '/',
		RES = '%',
		POW = '^',
		ASSIGN = '=',
		LP = '(',
		RP = ')'
	};

	/**
	* Класс для вычисления математических выражений
	* 
	* @author Bjarne Stroustrup :)
	* @author Gorobtsov Alexander
	* @author Getmanskiy Victor
	*/
	class Calculator
	{
		//int kpt,kpu,kmx,kse,kss,kspring;
		// Строка с разбираемым выражением
		const char* input_string;
		// Позиция в строке с выражением
		int cur_pos;
		// Таблица параметров
		std::vector<param> params_table;
		// Таблица параметров для сохранения
		std::stack<std::vector<param> > params_states;
		// Текущий TOKEN
		token_value current_token;
		// Текущее значение в выражении
		double number_value;
		// Текущее имя в выражении
		char name_string[256];

		// Поиск параметра в таблице по имени
		param* look(const char* name);

		// Устанавливает следующий TOKEN и возвращает его
		int get_token();
		// Вычислитель уровня 0
		double prim();
		// Вычислитель уровня 1
		double tterm();
		// Вычислитель уровня 2
		double term();
		// Вычислитель уровня 3
		double expr();

	public:

		Calculator();

		void Push();

		void Pop();

		sysparams sp;

		/** Вычислить выражение
		*
		* @param expression - вычисляемое выражение
		* @return вычисленное значение
		*/		
		double Eval(const char* expression);
	
		/** Вычислить выражение
		*
		* @param expression - вычисляемое выражение
		* @param value - вычисленное значение
		* @return true, если вычислено успешно
		*/
		bool TryEval(const char* expression, double& value);

		/** Вычислить выражение
		*
		* @param expression - вычисляемое выражение
		* @param value - вычисленное значение
		* @return true, если вычислено успешно
		*/
		double TryEval(const char* expression);

		/** Проверить вычислимость выражения
		* @return true, если вычислено успешно
		*/
		bool Check(const char* expression);
		
		// Добавить параметр-константу
		void AddParameter(const char* name, double value);
		
		// Добавить параметр-выражение
		void AddParameter(const char* name, const char* expression);

		// Пересчитать параметры
		void Update();
		// Пересчитать параметры до заданного включительно
		void UpdateTo(const char* name);
		// Пересчитать параметры после заданного включительно
		void UpdateAfter(const char* name);

		//Проверить наличие параметра
		int CheckParameter(const char* name, double* value = nullptr);
		
		//Вернуть последний добавленный в список параметр
		char* GetLastParameter();

		// Удалить добавленные после заданного параметры
		void RemoveAfter(const char* name);

		// Удалить параметр
		void RemoveParam(const char* name);

		// Добавить системный параметр
		double AddSystemParam(const char *s);

		// Назначить системные параметры
		//void SetSystemParams(int kpt,int kpu,int kmx,int kse,int kss = 0,int kspring = 0);

		//Проверить наличие параметра
		double GetParameter(const char* name);
		
		// Отладочный листинг таблицы параметров
		void Print();
	};
}


#endif // CALCULATOR_H
