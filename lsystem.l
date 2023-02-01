// Distributed System Transport Model

#include <lpfgall.h>
#include <math.h>

//defines pertinent to computation (INFTY, PI DT, etc)
#include "constants.h"

#include <list>
#include <iostream>

// include functions for Loading and Unloading, printing and A Hall's analytical solutions
#include "functions.h"

// defines pertinent to test case 

// constant source and sink
const float SOURCEj = 0.53;   	    // [g h-1] maximum initial influx from the source (for initialization of the system)
const float SINKj = 0.5;	    // [g h-1] maximun influx into the terminal sink
const float c0 = 0.4;               // [g cm-3] initial carbon concentration 
const float SOURCEcp = c0*c0/2.0;   //calculate carbon potential for c0

// defines pertinent to phloem resistance
const float R0 = 0.847E-3;	    // [MPa cm-4 h]
const float R = 0.024298245614035;  // [MPa cm3 g-1 K-1]
const float T = 293.15;             // [K]

/* -------------------------------------------------------------------------- */

// defines pertinent to water potential
const float Rx = 1.0063894024855E-6; 		// [cm-4 h] resistance per unit length of xylem
const float WP0 = 0.0;                          // [MPa]
const float MolalVolumeSucrose = 0.6301;        // [g-1 cm3]  
const double Jw = 0.0;                          // [cm3 h-1] 

// defines pertinent to system structure
const float INDE_LENGTH = 5.0;     	// [cm]
const int   START_SEGMENTS = 6;   	// number of segments to start on a shoot
const int   MAX_SEGMENTS = 60;    	// number of segments on a shoot

const float INT_SINKj = 0.005;    	// [g h-1] maximun possible influx into the internode sink
const float PRIORITY_INT = 0.01;  	// [g cm-3]
const float PRIORITY_SINK = 0.01; 	// [g cm-3]


/* -------------------------------------------------------------------------- */
// Declarations for carbon transport

// The structure representing a linearized organ/branch
struct LinearData
{
   float a; // intercept of the terminal element
   float b; // slope of the terminal element
};

// The structure characterizing a sink/source
struct SegmentData
{
  LinearData sink_source;   // linear approximation of the sink/source function {a,b} 
  float r;                  // internode resistance
  LinearData branch_flux;   // linear approximation of branch flux function 
  float j;                  // carbon flux into the sink/source part of the segment
  float cpt;                // carbon potential at the top node
  float cpb;                // carbon potential at the bottom node
  float err;                // current difference between iterations
  float C;                  // carbon in segment (structural or reserves)
  void Update_Branch_Flux();
  void Update_Carbon();
};
// Function incorporating flux from distal elements (Figure 4) adjusted for distributed sink
void SegmentData::Update_Branch_Flux()
{
  float denom = (1.0 + branch_flux.b*r +sink_source.b * r/2.);
  branch_flux.a = (branch_flux.a + sink_source.a + 
		   (sink_source.a * branch_flux.b - branch_flux.a * sink_source.b) * r/2.
		  )/denom;
  branch_flux.b = (branch_flux.b + sink_source.b)/denom;  
  return;
}


// Function to update the carbon in sink/source depending on carbon flow
void SegmentData::Update_Carbon()
{
  C += j * DT;
  if (C < 0.0)
    C = 0.0;
  return;
}

// Template for Segment initialisation to all zeros
SegmentData InitSegmentData = {{0.0, 0.0}, 0.0, {0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0}; 

/* -------------------------------------------------------------------------- */
// Declarations for water transport

// The structure used for water transport
struct HydraulicData
{
	double rs;      // hydraulic resistance
	double flow;    // water flow
	double pot;     // water potential 
	double Mcf;     // Munch counter flow
};

// template for hydraulic data initialisation
HydraulicData InitHydraulicData = {0.0,0.0,0.0,0.0}; 

/* -------------------------------------------------------------------------- */
// defines pertinent to simulation control

#define RUNS_TO_PRINT 1
const float END_SIM = DT * (MAX_SEGMENTS-START_SEGMENTS+5);


// Enumeration of computational phases
#define Develop         		1       // development phase
#define WaterFlowEst			21      // two phases for water flow
#define WaterPotential			22
#define Linearise       		31      // five phases for carbon dynamics 
#define UpdateBranchFlux		32
#define UpdateCarbonPotential           33
#define UpdateFlow      		34
#define Debug           		35
#define	Output                          41      // output phase

// Declarations of global variables 
int Phase;          // the current computation phase
float Serror;       // the error in computation for carbon dynamics, measuring convergence to the solution
int error_count;    // the number of times 'error' is above MAXERROR
float sim_time;     // the current simulation time (steps of DT)
float vis_time;     // the time to visualise the model

// Initialization of global variables
Start:
{
  Phase = Develop;
  Serror = 0.0;
  error_count = 0;
  sim_time = 0.0;
  vis_time = 0.0;

  // Set Output header information  
  FInit("output.csv",
        "Test parameter settings output\nInternodelen,Segments,c0, SOURCEcp , R0, Rx, Jw, SOURCEj, SINKj, INT_SINKj, PRIORITY_INT, PRIORITY_SINK\n");

  // output parameter values
  FPrintf("output.csv",
          "%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n\n",
       	  INDE_LENGTH, MAX_SEGMENTS, c0, SOURCEcp , R0, Rx, Jw,
          SOURCEj, SINKj, INT_SINKj, PRIORITY_INT, PRIORITY_SINK);

  // Header for Internode output parameter values
  FPrintf("output.csv",
          "Type,node,position,sim_time,sd_primary.r,sd_primary.j,sd_primary.cpt,h.pot,h.rs,h.flow,C\n"
         );

}

/* Processing at start of each derivation step */
StartEach:
{
  /* set rules to use for next phase */
  UseGroup(Phase);
  /* This switch statement controls derivation direction for the different phases of a model timestep */   
  /* specification of  phase computation is specified in the labelled groups follwiung the axiuom */
  switch (Phase)
  {
	// growth and development phase
    case Develop: 
      Forward();	            // left to right pass of the string
 
      if (sim_time < END_SIM)
      {	// update simulation time
        sim_time += DT;
      }
      else
      {
        Printf ("End of simulation!\n");
        Stop(); // stop simulation at end of run 
      }
      break;
	/* two phases for water flow */
    case WaterFlowEst:
		Backward();         // right to left pass
		break;
    case WaterPotential:            // left to right pass 
		Forward();
		break;
	/* five phases for carbon dynamics */
    case Linearise:                 // left to right pass
      Forward();
      break;
    case UpdateBranchFlux:          // right to left pass
      Backward();
      break;
    case UpdateCarbonPotential:     // left to right pass
      Forward();
      break;  
    case UpdateFlow:                // left to right pass
      Forward();
      Serror = 0;
      break;
    case Debug:                     // prints out info of all sinks/sources with error > MAX_ERROR
      Printf ("Warning: error > MAXERROR (%g > %g) after %d iterations\n", Serror, MAXERROR, error_count);
      Printf ("Dumping 'S' modules at sim_time %g:\n", sim_time);
      Forward();
      break;
    case Output:                    /* one phase for output */
	    Printf ("Start of Output Phase at sim_time %g:\n", sim_time);
		Forward();
		break;
  }
}

/* Processing at end of each derivation step  */
EndEach:
{   
  // this switch statement determines next phase type
  switch (Phase)
  { 
    // growth and development phase
    case Develop:
      if (sim_time <= END_SIM)
      {
          Phase = WaterFlowEst;	    // switch to hydraulic processing
          // Set to display frame if appropriate
          if (sim_time >= vis_time)
          {
            DisplayFrame();
            vis_time += VIS_STEP;
          }
          error_count = 0;
      }
      break;
	/* two phases for water flow */
    case WaterFlowEst:
      Phase = WaterPotential;
      break;
    case WaterPotential:
      Phase = Linearise;	   // switch to carbon processing 
      break;
	/* five phases for carbon dynamics */
    case Linearise:
      Phase = UpdateBranchFlux;
      break;
    case UpdateBranchFlux:
      Phase = UpdateCarbonPotential;
      break;
    case UpdateCarbonPotential:
      Phase = UpdateFlow;
      break;
    case UpdateFlow:
      if (Serror < MAXERROR)      // successful flow update processing
      {
        if (sim_time >= END_SIM-DT*RUNS_TO_PRINT) {  // change to output phase at end of run
          Phase = Output;
        } else { // continue to next step
          Phase = Develop;
        }
        error_count = 0;
      }
      else  // further processing necessary
      {
        if (error_count > 1000) // too many passes? 
        {	// report problem
          Phase = Debug;
        }
        else
        {	// iterate process
          ++error_count;
          Phase = Linearise;
        }
      }
      break;
    case Debug:                 /* stop after debug */
      Stop();
      break;
    case Output:                /* one phase for output */
	  Printf ("End of Output Phase at time %g:\n", sim_time);
	  Phase = Develop;
	  break;
  }
}

/* -------------------------------------------------------------------------- */
// module type declarations
module S(SegmentData);          // a sink/source component
module Hd(HydraulicData);       // a hydraulics component
module Internode(int,float);    // internode with attributes: node,length 
module Sink(int);            	// sink with attributes: node
module TestCase();              // the initial setup of structure
/* -------------------------------------------------------------------------- */

derivation length: 999999;

Axiom: TestCase();

decomposition:

TestCase():
{
  // each segment of shoot is modelled as one internode
  SegmentData sd_internode_primary = InitSegmentData;
  SegmentData sd_sink = InitSegmentData;
  SegmentData sd_source = InitSegmentData;
 
  HydraulicData h_internode = InitHydraulicData;
  HydraulicData h_sink = InitHydraulicData; 

  h_internode.rs = Rx * INDE_LENGTH; 	// hydraulic resistance
  h_internode.pot = -Jw*Rx*INDE_LENGTH;	// initial hydraulic potential

  // set the resistance along transport pathway (in the first S module associated with Internode)
  const float Rs = ( R0 * PhloemResistanceRatio (c0)/(PressureRatio (c0)*R*T));   // resistance per unit length (initial case)
  sd_internode_primary.r = Rs*INDE_LENGTH;  // based on Alistair's parameters  
 
    // source at base is defined as a boundary condition in UpdateCarbonPotential code
    /* produce segments */
    for (int node=1; node<=START_SEGMENTS; node++) { 
		// initial conditions for cp top and bottom
		sd_internode_primary.cpt = SOURCEcp - node * INDE_LENGTH * (SOURCEj - INT_SINKj * node * INDE_LENGTH/2) * Rs;    
 		sd_internode_primary.cpb = SOURCEcp - (node-1) * INDE_LENGTH * (SOURCEj - INT_SINKj * node * INDE_LENGTH/2) * Rs;

 		nproduce 
          	  Internode(node,INDE_LENGTH)
              S(sd_internode_primary)
 	    	  Hd(h_internode)
		      // to ensure consistency for rule matching, place SB() EB() here, where a leaf would be.
              SB() /* Leaf */ EB()
         ;
    }
    /* produce sink */
    sd_sink.j = SINKj;  // initialise for test case 1 
    sd_sink.cpt = SOURCEcp - MAX_SEGMENTS * INDE_LENGTH * SOURCEj * Rs; 
    sd_sink.cpb = SOURCEcp - MAX_SEGMENTS * INDE_LENGTH * SOURCEj * Rs;
    h_sink.flow = Jw;
    produce
            Sink(START_SEGMENTS+1)
            S(sd_sink)
	        Hd(h_sink)
	       // to ensure consistency for rule matching, place SB() EB() here, where a leaf would have been.
            SB() /* Leaf */ EB()
          ;
}

production:

/* Development Rules -------------------------------------------------------- */

group Develop: // update of parameters


Internode(node,length) S(sd_primary)  Hd(h): 
{
  // accumulate carbon for primary growth
  sd_primary.Update_Carbon();
  double carbon = sqrt(sd_primary.cpb + sd_primary.cpt);
// set resistance along transport pathway based on local C concentration (derived from average carbon potential)
  sd_primary.r = R0 * PhloemResistanceRatio(carbon) / (PressureRatio(carbon)*R*T)*INDE_LENGTH;


 // Effect of phloem flux on hydraulic flow: Munch counter flow (following Equation 5) 
  double branchj = sd_primary.branch_flux.a + sd_primary.branch_flux.b * sd_primary.cpb;
  h.Mcf = MolalVolumeSucrose * branchj - branchj / carbon;
 		
  produce Internode(node,length) S(sd_primary) Hd(h);
}

/* Update sink, if necessary, create new segment before sink with same content as previous segment */
Internode(inode,ilength) S(sd_internode_primary) Hd(h_internode) << Sink(node) S(sd_sink)  Hd(h):
{
  // update sink size
  sd_sink.Update_Carbon();

  if (inode<MAX_SEGMENTS) {// produce new internodes every develop step, until maximum number reached
    nproduce 
          Internode(inode+1,ilength)
          S(sd_internode_primary)
 	  Hd(h_internode)
          SB() /* Leaf */ EB() // to ensure consistency for rule matching, place SB() EB() here, where a leaf would have been.
         ;
    node = node+1;
    Printf("new node: %d\n", inode+1);
  }

  produce Sink(node) S(sd_sink)  Hd(h);
}

S(sd):
{
  // if r == INFTY this is an inert S between branches
  if (sd.r != 0.0 || sd.sink_source.a != 0.0 || sd.sink_source.b != 0.0)
    Printf ("PANIC! Orphan S() module found. Check 'Develop' phase.\n");  
  produce S(sd);
}


/* -------------------------------------------------------------------------- */


/* Carbon Dynamics Phaes defined as groups ---------------------------------- */

group Linearise: // linearise equations for sinks/sources 


Internode(node,length) S(sd_primary) :
{
  //carbon flux function is set
  sd_primary.sink_source.a = INT_SINKj *
        (SinkUnloading(sd_primary.cpt, PRIORITY_INT)
	-sd_primary.cpt	* dSinkUnloading_dcp (sd_primary.cpt, PRIORITY_INT));

  sd_primary.sink_source.b = INT_SINKj *
	dSinkUnloading_dcp (sd_primary.cpt, PRIORITY_INT);

  produce Internode(node,length) S(sd_primary) ;
}


Sink(node) S(sd_sink) :
{
  sd_sink.sink_source.a = SINKj *
        (SinkUnloading(sd_sink.cpt, PRIORITY_SINK)
	-sd_sink.cpt * dSinkUnloading_dcp (sd_sink.cpt, PRIORITY_SINK));

  sd_sink.sink_source.b = SINKj *
	dSinkUnloading_dcp (sd_sink.cpt, PRIORITY_SINK);
  
  produce Sink(node) S(sd_sink);

}

S(sd):
{
  // if r, a and b == 0 this is an inert S between branches
  if (sd.r != 0.0 || sd.sink_source.a != 0.0 || sd.sink_source.b != 0.0)
    Printf ("PANIC! Orphan S() module found. Check 'UpdateDemand' phase.\n");
  produce S(sd);
}

/* -------------------------------------------------------------------------- */
// calculate equivalent sink/source function of the  branch starting at the base

group UpdateBranchFlux:	

consider: S;

S(sd) >> SB() S(sdr2) EB() S(sdr1):
{ 

  sd.branch_flux.a = sdr1.branch_flux.a+sdr2.branch_flux.a;
  sd.branch_flux.b = sdr1.branch_flux.b+sdr2.branch_flux.b;

  sd.Update_Branch_Flux();
  produce S(sd);
}

S(sd) >> SB() S(sdr) EB():
{
  sd.branch_flux.a = sdr.branch_flux.a;
  sd.branch_flux.b = sdr.branch_flux.b;

  sd.Update_Branch_Flux();
  produce S(sd);
}


S(sd) >> S(sdr):
{
  sd.branch_flux.a = sdr.branch_flux.a;
  sd.branch_flux.b = sdr.branch_flux.b;

  sd.Update_Branch_Flux();
  produce S(sd);
}


S(sd):
{  
  sd.branch_flux.a=0.0;
  sd.branch_flux.b=0.0;

  sd.Update_Branch_Flux(); 
  produce S(sd);
}

/* -------------------------------------------------------------------------- */

group UpdateCarbonPotential: // find the concentration of carbon along transport pathway

consider: S;

// calculate the carbon potentials; estimate the error
S(sdl) << S(sd):
{
  sd.cpb = sdl.cpt;
  sd.cpt = (sd.cpb*(1.0 - sd.branch_flux.b*sd.r) - sd.branch_flux.a*sd.r + sd.sink_source.a*sd.r/2.)/
           (1- sd.sink_source.b*sd.r/2.); 
  produce S(sd);
}

S(sd):
{  // Boundary condition sets carbon potential at the bottom equal to SOURCEcp
  sd.cpb = SOURCEcp;
   // calculation for carbon potential at the top
  sd.cpt = (sd.cpb*(1.0 - sd.branch_flux.b*sd.r) - sd.branch_flux.a*sd.r + sd.sink_source.a*sd.r/2.)/
           (1- sd.sink_source.b*sd.r/2.); 
  produce S(sd);
}



/* -------------------------------------------------------------------------- */

group UpdateFlow: // update the flow of carbon into each sink/source

S(sd) :
{ 
  extern float Serror;
  float new_j, divisor;

  if (sd.cpt >= 0.0)
  {
    // Calculate new flux into sink_source
    new_j = sd.sink_source.a + sd.sink_source.b * sd.cpt; 

    // Calculate error and save for display purposes
    sd.err = new_j - sd.j;

    // Calculate cumulative error to determine whether it is possible to stop
    Serror += fabs(sd.err);

    // Update flux to be consistent with present values of cpt
    // check for minimum current
    if (fabs(new_j) < EPSILON)
      sd.j = 0.0;
    else
      sd.j = new_j;
  }
  else // cpt < 0! Reset the Newton-Raphson method, with cpt as zero
  {
    sd.cpt = 0.0;
    sd.err = 1.0;
    Serror += 1.0;
  }

  produce S(sd);
}

/* -------------------------------------------------------------------------- */

group Debug: // if error doesn't fall below MAXERROR after many iterations, output debug info

S(sd):
{
  if (sd.err > MAXERROR)
    Printf ("C:%g j:%g cpt:%g err:%g\n", sd.C, sd.j, sd.cpt, sd.err);

  produce S(sd);
}


/* --------------------------------------------------------------------------- */
/* Accumulate water flows from right to left                                   */

group WaterFlowEst:

consider: Hd;

Hd(h) >> SB() EB() SB() Hd(hr):
{
	h.flow = hr.flow;
	
	produce Hd(h);
}


Hd(h) >> SB() Hd(hr2) EB() Hd(hr1):
{
	h.flow = hr1.flow + hr2.flow;
	
	produce Hd(h);
}


Hd(h) >> SB() Hd(hr) EB():
{
	h.flow = hr.flow;
	
	produce Hd(h);
}


Hd(h) >> Hd(hr):
{
	h.flow = hr.flow;
	
	produce Hd(h);
}


Hd(h):
{
	produce Hd(h);
}

/* ----------------------------------------------------------------------------- */
/* Calculate water potential from left to right incorporating Munch counter flow */


group WaterPotential:

Hd(hl) << Hd(h):
{
	h.pot = hl.pot - (h.flow + h.Mcf) * h.rs;
	
	produce Hd(h);
}

Hd(h):
{
	h.pot = -(h.flow + h.Mcf) * h.rs;

	produce Hd(h);
}

/* -------------------------------------------------------------------- */
/* at every OUTPUT_STEP time step, write model output to a file         */

group Output: 


Internode(node,length) S(sd_primary) Hd(h) : 
{
  FPrintf("output.csv","Int,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", 
           node, node*length, sim_time, sd_primary.r, sd_primary.j, sd_primary.cpt,
	   h.pot, h.rs, h.flow,
	   sqrt(2.0*sd_primary.cpt)
	 );
  
  produce Internode(node,length) S(sd_primary)  Hd(h) ; 
}

Sink(node) S(sd_sink) Hd(h):
{
  
  FPrintf("output.csv","Snk,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", 
          node, (node-1)*INDE_LENGTH, sim_time, sd_sink.r, sd_sink.j, sd_sink.cpt,
          h.pot, h.rs, h.flow,
          sqrt(2.0*sd_sink.cpt)
         );
  produce Sink(node) S(sd_sink) Hd(h);

}

/* -------------------------------------------------------------------------- */
/* Default Group ------------------------------------------------------------ */

group 0:

/* Visualisation ------------------------------------------------------------ */

interpretation:
maximum depth: 2;


Internode(node,length)  S(sd_primary)  Hd(h): 
{
  int colour = 15;  
 
  produce SetColor(colour) SetWidth(1) F(length);
}

Sink(node):
{
  float radius=.5;
  produce SB() 
          SetColor(13) f(radius) Sphere(radius)
          EB();
}


/* -----lsystem.l--------------------------------------------------------------------- */
