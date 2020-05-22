#ifndef _FLAG_H_
#define _FLAG_H_

#include "utils.h"
#include <vector>
void updateFlagField(Fields &fields, std::vector<Cell> &filledCells, std::vector<Cell> &emptiedCells, int * length, int n_threads);

#endif
