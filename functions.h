/* This file contains print functions, sink/source functions, phloem resistance and osmotic pressure functions */
/* If you make any changes, don't forget to save it before running the model. */

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H
#include <cstdarg>
// File print functions
void FInit (const char* filename, const char * format, ... )
{
	FILE* stream = fopen(filename,"w+");
	
	va_list args;
	va_start (args, format);
	vfprintf (stream, format, args);
	va_end (args);
	fclose(stream);
}


void FPrintf (const char* filename, const char * format, ... )
{
	FILE* stream = fopen(filename,"a+");
	
	va_list args;
	va_start (args, format);
	vfprintf (stream, format, args);
	va_end (args);
	fclose(stream);
}


/* -------------------------------------------------------------------------- */
// Common sink/source functions
/* -------------------------------------------------------------------------- */

float SinkUnloading (float cp, float q) 
// Michaelis-Menton function expressed in cp, defines carbon unloading
// cp is carbon potential cp=c^2/2, (g^2 C / cm^6)
// q is priority parameter (g C / cm^3)
{
  float c = sqrt(2.0*cp); // carbon concentration at the top of the node
  return ((c / (q + c)));
}

float dSinkUnloading_dcp (float cp, float q) 
// defines the derivative  of carbon unloading function with respect to 'cp',  
//              df/dcp=(df/dcp)*(dcp/dv)=(df/dcp)*(1/cp)
// cp is carbon potential cp=c^2/2, (g^2 C / cm^6)
// q is priority parameter (g C / cm^3)
{
  float c = sqrt(2.0*cp);
  return (q / ((q + c)*(q + c) * (c + 0.001))); // added 0.001 to avoid dividing by 0.
}

float SourceLoading (float cp, float q) 
// cp is carbon potential cp=c^2/2 (g^2 C / cm^6)
// q is loading rate parameter
{
  float c = sqrt(2.0*cp);
  // source loading = -1 + (c / (q + c))
  return (-1.0 + c / (q + c));
}

float dSourceLoading_dcp (float cp, float q) 
// defines the derivative of carbon loading function with respect to 'cp' 
{
  float c = sqrt(2.0*cp);
  return (q / ((q + c) * (q + c) * (c + 0.001)));
}


/* -------------------------------------------------------------------------- */
 /* Functions relevant to phloem resistance and osmotic pressure */
/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------*/
// phloem resistance ratio R/R0     
float PhloemResistanceRatio (float cab) 
{
  float  k1 = 0.19961;
  float  k2 = 1.00954;
  float  k3 = -2.2249;
  float  k4 = 1.2776;
 
  return (1.0+k1*cab+k2*cab*cab)/(1.0+k3*cab+k4*cab*cab);			 
} 

/* -------------------------------------------------------------------------- */
//Osmotic pressure ratio P/RT

float PressureRatio (float cab) 
{
  float  a = 1.03913;
  float  b = 0.86536;
  float  c = 5.2465;
  
  return (a + b*cab + c*cab*cab);			 
} 

#endif
