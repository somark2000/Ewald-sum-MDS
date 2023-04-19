#pragma once

#include "datatypes.h"

#define ComPair(i1, j1, i2, j2) \
   ((i1 == i2 && j1 == j2) || (i1 == j2 && j1 == i2))


int GetPairNum0(Atom atomsf[], int natm, Bond bndsf[], int &nbnd);
void PairList0(Atom atoms[], int natm, Bond bnds[], int nbnd, Pair pairs[], int &npair);