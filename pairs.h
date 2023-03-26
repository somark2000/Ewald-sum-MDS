#pragma once

#include "datatypes.h"

#define ComPair(il, jl, i2, j2) ((il == i2 && jl == j2) || (il == j2 && jl == i2))

int GetPairNum0(Atom atomsf[], int natm, Bond bndsf[], int nbnd);
void PairList0(Atom atoms[], int natm, Bond bnds[], int nbnd, Pair pairs[], int& npair);