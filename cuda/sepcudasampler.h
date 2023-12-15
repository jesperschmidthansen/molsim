
#ifndef __SEPCUDASAMPLER_H__
#define __SEPCUDASAMPLER_H__

#include "sepcudadefs.h"
#include "sepcudamisc.h"

typedef struct {
	
	double **dacf;
	double **tmacf;
	double **stress;
	
	double **mcoskrArray;
	double **msinkrArray;

	double **vcoskrArray;
	double **vsinkrArray;
	
	double **stressa, **stressb;
	
	double *wavevector;
	
	unsigned int lvec, nwaves;
	unsigned int index, nsample;
	
	double dtsample;
} sepcugh;

typedef struct {
	
	double **stress;
	double **stressax, **stressbx;
	double **stressay, **stressby;
	double **stressaz, **stressbz;
	unsigned int stresslvec, stressindex, stressnsample;

	double **dipole;
	double **dipoleax, **dipolebx;
	double **dipoleay, **dipoleby;
	double **dipoleaz, **dipolebz;
	unsigned int dipolelvec, dipoleindex, dipolensample;

	unsigned nwaves;
	
	double *wavevector;
	double dtsample;

} sepcumgh;

typedef struct {
	double **corr;
	double **dipoleax, **dipolebx;
	double **dipoleay, **dipoleby;
	double **dipoleaz, **dipolebz;
	unsigned int lvec, index, nsample;

	unsigned nwaves;
	double *wavevector;
	double dtsample;
} sepcusampler_dipole;

typedef struct {
	double **corr;
	double **stressa, **stressb;

	unsigned int lvec, index, nsample;

	unsigned nwaves;
	double *wavevector;
	double dtsample;
} sepcusampler_stress;

// Aux
double** sep_cuda_matrix(size_t nrow, size_t ncol);
void sep_cuda_free_matrix(double **ptr, size_t nrow);

// gh samlper
sepcugh* sep_cuda_sample_gh_init(sepcusys *sysptr, int lvec, unsigned nk, double dtsample);
void sep_cuda_sample_gh(sepcugh *sampleptr, sepcupart *pptr, sepcusys *sptr);
void sep_cuda_sample_gh_free(sepcugh *ptr);

// mgh sampler - depreciated
void sep_cuda_print_current_corr(sepcusys *sptr, sepcumgh *sampler, const char quantity, const char *filename);
sepcumgh* sep_cuda_sample_mgh_init(sepcusys *sysptr, int lvec[2], unsigned nk, double dtsample);
void sep_cuda_sample_mgh_free(sepcumgh *ptr);
void sep_cuda_sample_mgh(sepcumgh *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr);

// Polarization correlations
sepcusampler_dipole* sep_cuda_sample_dipole_init(sepcusys *sptr, int lvec, unsigned nk, double dtsample);
void sep_cuda_sample_dipole_free(sepcusampler_dipole *ptr);
void sep_cuda_sample_dipole(sepcusampler_dipole *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr);

// Stress correlations
sepcusampler_stress* sep_cuda_sample_stress_init(sepcusys *sysptr, int lvec, unsigned nk, double dtsample);
void sep_cuda_sample_stress_free(sepcusampler_stress *ptr); 
void sep_cuda_sample_stress(sepcusampler_stress *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr);


#endif
