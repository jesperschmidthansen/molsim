#include "sep.h"

typedef struct{
  unsigned int npart;
  
  char type;
  unsigned int npart_type;

  double tnow, dt; 

  double *m, *z;
  double **x, **v, **s, **w, **inertia;
} data;

void allocate_memory(data *dptr);
void free_memory(data *dptr);

int file_exists(const char file[]);
void read_entry(data *dptr, const char file[]);
void get_time(data *dptr, const char file[]);
void get_num_type(data *dptr, const char type, const char file[]);
