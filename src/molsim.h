
#ifndef __MOLSIM_H__

#include "mex.h"
#include "sep.h"
#include "task.h"
#include <string.h>

#define MAXNUMBTASK 12

// Globals static variables *sigh*
static sepatom *atoms;
static sepsys sys;
static sepret ret;
static sepmol *mols;
static sepsampler sampler;
static taskmanager *tasks;

static int natoms;
static long unsigned int iterationNumber = 0;
static double alpha[3] = {0.1};    
static int exclusionflag = SEP_ALL;
double SPRING_X0=0.0; 

static bool initflag = false;
static bool tempflag = false;
static bool initmol = false;
static bool initsampler = false;
static bool inittasks = false;

static double lbox[3], dt=0.005, maxcutoff=2.5, temperature=1.0,
  compressionfactor = 0.9995, taufactor=0.01; 

static unsigned int ntasks = 0;

// Hard-coded hash values for switch - *I* cannot "optimize" further
// Hash value is simply the string (lower case) character sum
enum {
  RESET=547, CALCFORCE=930, INTEGRATE=963,
  THERMOSTATE=1200, SAMPLE=642, ADD=297,
  GET=320, PRINT=557, SAVE=431,
  TASK=435, COMPRESS=876, CLEAR=519,
  SET=332, HELLO=532, LOAD=416,
  HASHVALUE=961, BAROSTATE=965
};

// Wrapper functions for actions
void action_reset(int nrhs);
void action_set(int nrhs, const mxArray **prhs);
void action_load(int nrhs, const mxArray **prhs);
void action_calcforce(int nrhs, const mxArray **prhs);
void action_thermostate(int nrhs, const mxArray **prhs);
void action_barostate(int nrhs, const mxArray **prhs);
void action_integrate(int nrhs, const mxArray **prhs);
void action_save(int nrhs, const mxArray **prhs);
void action_print(void);
void action_get(mxArray **plhs, int nrhs, const mxArray **prhs);
void action_sample(int nrhs, const mxArray **prhs);
void action_task(int nrhs, const mxArray **prhs);
void action_compress(int nrhs, const mxArray **prhs);
void action_clear(int nrhs, const mxArray **prhs);
void action_add(int nrhs, const mxArray **prhs);
void action_hash(int nrhs, const mxArray **prhs);

// Local helper functions
double spring_x0(double r2, char opt);
void inputerror(const char* funstr);
unsigned hashfun(const char *key);

// Extern
void sep_lattice(int, char **);
void sep_sfg(int, char **);


#endif
