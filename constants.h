#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// defines pertinent to computation
const double INFTY     = 1e10;      // very large number
const double EPSILON   = 1.0/INFTY; // very small number different from zero
const double DT        = 0.04166666666666666666667;       // time step (days), 1 hour
//const double DT        = 0.1;       // time step (days)
const double MAXERROR  = 1e-4;      // error threshold for evaluating carbon flow
const double MAXERRORH = 1e-2;      // error threshold for evaluating carbon flow
const double ENV_STEP  = 1.0;       // when to compute light distribution (days) (keep as 1 for now!!!)
const double VIS_STEP  = 1.0;       // when to visualise model (days)
const double OUT_STEP  = 0.04166666666666666666667;      // when to write model output to a file
const double PI        = 3.14159265358979323846264; // PI

#endif
