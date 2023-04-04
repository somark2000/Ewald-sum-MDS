#pragma once

#include <omp.h>
#include "basetype.h"
#include "Constants.h"
#include "datatypes.h"
#include "Memalloc.h"

void EwaldK0(Atom atoms[], int natm, VecR3 Box, real alpha, real& Eele);
