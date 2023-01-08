#pragma once

#include "basetype.h"

//basic mathematic constants
const real pi = 3.141592653589793;
const real rad = pi / 180;
const real deg = 180 / pi;

//basic physics constants 
const real amu = 1.660539040e-27;
const real NA = 6.022140857e23;
const real kB = 1.38064852e-23;
const real e = 1.6021766208e-19;
const real eps0 = 8.854187817e-12;
const real c = 299792458e0;
const real kcal = 4184e0;

//basic conversion constants
const real A = 1e-10;
const real nm = 1e-9;
const real ps = 1e-12;
const real eV = e;
const real bar = 1e5;
const real facc = kcal / (NA * amu) * (ps * A)*(ps * A);
const real fkin = 1e0 / facc;
const real fcoul = e * e * NA / (4e0 * pi * eps0 * A * kcal);