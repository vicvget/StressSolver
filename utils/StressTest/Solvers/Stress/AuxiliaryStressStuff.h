#ifndef AUXILIARY_STRESS_STUFF_H

#define AUXILIARY_STRESS_STUFF_H


#include <bitset>
#include <vector>


// количество степеней свободы у одного элемента, на которые разбивается тело
// при расчете напряженно-деформированного состояния
const size_t FreedomsCount = 6;


// набор флагов, закреплена ли соответствующая степень свободы
// в некотором элементе сеточного представления тела
using IsSealedFlags = std::bitset<FreedomsCount>;


// список наборов флагов закрепленных степеней свободы для элементов сеточного представления тела
using IsSealedFlagsList = std::vector<IsSealedFlags>;


#endif // AUXILIARY_STRESS_STUFF_H