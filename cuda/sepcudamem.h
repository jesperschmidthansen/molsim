
#ifndef __SEPCUDAMEM_H__
#define __SEPCUDAMEM_H__


#include "sepcudadefs.h"
#include "sepcudamisc.h"

sepcupart* sep_cuda_allocate_memory(unsigned npart);
sepcupart* sep_cuda_load_xyz(const char *xyzfile);
void sep_cuda_free_memory(sepcupart *ptr);
sepcusys *sep_cuda_sys_setup(sepcupart *pptr);


#endif
