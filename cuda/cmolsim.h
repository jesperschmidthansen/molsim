/**************************************
 *
 * Header for the cmolsim
 *
 **************************************/

void load_xyz(const char xyzfile[]);
void load_top(const char topfile[]);
void save_xyz(const char filename[]);

void free_memory(void);

void reset_iteration(void);

void force_lj(const char *types, float *ljparam);
void force_coulomb(float cf);
void force_bond(int type, float lbond, float ks);
void force_angle(int type, float angle, float kangle);
void force_torsion(int type, float *params);

void integrate_leapfrog(void);

void thermostat_nh(const char *type, float temp0, float thermostatmass);
void reset_momentum(int resetfreq);

void get_pressure(double *press);
void get_energies(double *energies);
void get_positions(double *positions);
void get_velocities(double *velocities);
void get_masses(double *masses);
void get_charges(double *charges);
void get_types(char *types);
unsigned get_npart(void);

void set_exlusion_molecule(const char rule[]);
void set_timestep(float dt);
void set_maximum_cf(float cf);

