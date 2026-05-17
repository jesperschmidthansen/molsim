 <html>
<body>

UNDER DEVELOPMENT 

This text is not meant to introduce molecular dynamics; such introductions can be found in many standard books. 
In brief, the basic idea is to solve the classical 
equation of motion of an ensemble of interacting particles. In the simplest form this means solving (numerically) Newton's second law

$$
   \frac{\mathrm{d}\mathbf{r}_i}{\mathrm{d} t} = \mathbf{v}_i \ , \ \
   \frac{\mathrm{d}\mathbf{p}_i}{\mathrm{d} t} = \mathbf{f}_i \ ,
$$

where $\mathbf{r}_i, \mathbf{v}_i, \mathbf{p}_i$ and $\mathbf{f}_i$ are the particle position, velocity, momentum and force acting on the particle. In a
standard simulation we solve this set of differential equations by (i) evaluating the forces acting on the particles, and (ii) from this integrate forward in time. The following pseudo code lists the basic idea

<pre><code>
Set initial configuration: positions r and momenta p
 
do (as many times as we want)
   f <- calcforce(r)
   r, p <- integrate(f,p)
done
</pre></code>
<p> The loop is here denoted the main MD-loop.</p>

Below you can see how this is implemented in molsim
<pre><code>
# Instance of molsim object
sim = molsim();

# Set number of particles 10x10x10 and simulation box lengths to 10.557 in all three directions. 
# Temperature is initially set to 1.0
sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);

# Main MD loop
for n=1:1e4
   # Calculates the forces acting on particles of type A - default type
	 sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   
   # Integrate forward in time using the leap-frog algorithm
	 sim.leapfrog();
end
</code></pre>

Please see the example folder, where different and more complicated simulations are shown.

<h3>Force field model</h3>
The force is given by the gradient of the potential energy function $U$ by $\mathbf{f} = - \nabla U$. 
Currently, the energy function in molsim is given by 

$$
 U(\mathbf{r}_i, r_{ij}, \ldots) =  U_\mathrm{lattice} + U_\mathrm{vWaals} + U_{\mathrm{coulomb}} + U_\mathrm{bonds} + U_\mathrm{angles} +  U_\mathrm{torsion}
$$

The table shows the terms 
<table>
 <tr> <td> Potential function </td> <td> User supplied parameters </td> <td> Method </td> </tr>
 <tr>
  <td> $U_\mathrm{lattice} =  \frac{1}{2}\sum_\mathrm{sites}k_0 (\mathbf{r}_i - \mathbf{r}_\text{lattice})^2$  </td>
 <td> $k_0$</td>
 <td> atoms.tether(atom type, $k_0$) </td>
 </tr>
 <tr>
  <td>$U_\mathrm{vWaals} =  4\sum_{i,j \, \mathrm{pairs}} \epsilon\left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - a_w \left(\frac{\sigma}{r_{ij}}\right)^{6}\right]$ </td>
  <td> $r_\text{cutoff}$, $\epsilon$, $\sigma$, $a_w$</td>
  <td> lennardjones(atoms types, [ $r_\text{cutoff}$, $\epsilon$, $\sigma$, $a_w$ ]) </td>
 </tr> 
<tr> 
 <td> $U_{\mathrm{coulomb}} = \sum_{i,j \, \mathrm{pairs}}\frac{q_iq_j}{r_{ij}}$</td> 
<td> $r_\text{cutoff}$ </td>
<td> sfcoulomb($r_\text{cutoff}$) </td>
</tr>
<tr><td> $U_{\mathrm{bonds}} =\frac{1}{2} \sum_{\mathrm{bonds}} k_{s}(r_{ij} - l_0)^2$ </td>
<td> $k_s$, $l_0$ </td>
<td> harmonicbond() (params set before call) </td>
</tr>
<tr> 
 <td> $U_{\mathrm{angles}}=\frac{1}{2}\sum_{\mathrm{angles}} k_{\theta} (\cos(\theta) - \cos(\theta_0))^2$ </td>
<td> $k_\theta$, $\theta_0$ </td>
<td> cossqangle() (parameters set before call)</td>
</tr>
<tr> 
 <td> $U_\mathrm{torsion}=\sum_{\mathrm{angles}} \sum_{n=0}^5 c_n \cos^n(\pi-\phi)$ </td>
 <td> $c_n$ </td>
 <td>ryckbell() (parameters set before call)</td>
</tr>
</table>

Notice that the different terms can be mapped to and from other force fields or have approximated similar behavior around the minimum energies. 

<h3>Integrator and thermostats</h3>
Currently, molsim only includes the leap-frog integrator. The call is simply
<pre><code>
molsim.leapfrog();
</code></pre>
The integrator time step is changed by accessing the integrator member directly; typically before the main MD-loop 
<pre><code>
sim = molsim();
sim.integrator.dt = 5e-4; # Default is 5e-3
</code></pre>

molsim includes two different thermostats. The Nose-Hoover thermostat and a simple relaxation-type
thermostat. The thermostat type is set before the main MD-loop
<pre></code>
molsim.setthermostat(type of thermostat, atom types, target temperature, relaxation parameter);
molsim.setthermostat(type of thermostat, target temperature, relaxation parameter);
</code></pre>
The thermostat type can be "nh" or "relax". If atom type is not specified, the thermostat is applied
to all atoms. See method help text for relaxation parameter. To apply the thermostat you call 
<pre><code>
molsim.applythermostat();
</code></pre>
in the main MD-loop. 

<h3>Setting and saving the configuration</h3>
There are two ways to set the (initial) configuration
<pre><code>
molsim.setconf([Nx, Ny, Nz], [Lx, Ly, Lz], temperature);
molsim.setconf(filename);
</code></pre>
1: In the first way, the user specifies the number of atoms and simulation box length in each of
the three dimensions, as well as provide an initial kinetic temperature which defines the initial
atom velocities. Mass, type, and charge are set to 1, 'A', and 0, respectively; these can be set
after the call to setconf.

2: The user can set the configuration from a file; this file must be in either xyz-format or
mat-format. If the mat-format is used, it is strongly recommended that the configuration has been
saved using molsim's save functionality.

To save a configuration you can use the save method
<pre></code>
molsim.save(filename);
</code></pre> 
where filename as the extension .xyz or .mat. The mat-format is recommended as this saves additional
informations that are useful for data analysis. The xyz-format
enables the user to visualize the system with external tools like VMD and Ovito.

You can also autosave;  the autosaver saves simulation snaps shot of time, configuration space, 
force and simulation box crossing. This is useful for post simulation data analysis. The autosaver outputs 
a compressed file 'molsim.zip' containing data files  'molsim-%06.mat" with variables 
'time', 'r', 'v', 'f', 'bxcrs'.  If a file 'molsim.zip' exists in the current directory it will be deleted.
The autosaver will slow down the execution time. The autosaver is set before the main MD-loop

<pre><code>
sim = molsim();
sim.setautosave(100); # Save every 100 time steps
</code></pre>


<h3>Using molconf</h3>


<h3>Why MEX?</h3>
<p>GNU Octave offers a great C++ interface with the dynamically linked functions (DLDs). However, my experience is that the pure C MEX interface to produces faster running binaries. This is perhaps due to DLD's call-by-value interface giving an additional copying overhead.    
</p>

<p>
Test: The functions below shows a DLD and a MEX version of a function that calculates the sum of an array and updates the array with a number; this is a relevant task in molecular dynamics. 
<table>
 <tr> <td> DLD msum_oct.cpp </td><td>MEX msum_mex.c</td></tr>
<tr>
 <td>
  
 <pre><code>
#include &lt;octave/oct.h&gt;

DEFUN_DLD(msum_oct, args, ,""){
   octave_value_list retval;
   Matrix A(args(0).array_value());
   int nrows = A.dim1();
   int ncols = A.dim2();

   double *Aptr = A.fortran_vec();

   double sum=0.0f;
   for (int n=0; n&lt;nrows; n++) {
     for (int m=0; m&lt;ncols; m++) {
         int idx = m*nrows + n;
         sum += Aptr[idx];
         Aptr[idx] += 1.0;
      }
   }

   retval.append(sum);
   retval.append(A);
   return retval;
}
</code></pre>
</td>

<td>
 
<pre><code>
 #include "mex.h"
 
 void mexFunction(int nlhs, mxArray *plhs[], 
              int nrhs, const mxArray *prhs[]) {
 
     double *A = mxGetPr(prhs[0]);
     int nrows = mxGetM(prhs[0]);
     int ncols = mxGetN(prhs[0]);
 
     double sum=0.0f;
     for (int n=0; n&lt;nrows; n++) {
        for (int m=0; m&lt;ncols; m++) {
           int idx = m*nrows + n;
           sum += A[idx];
           A[idx] += 1.0;
        }
      }
 
      plhs[0] = mxCreateDoubleScalar(sum);
 }

 
 
</code></pre>
</td>
</table>
</td>
</p>

<p>The functions are compiled with or without <code>-Ofast</code> flag. Timing is done by 
<pre><code>
>> A=randn(1000, 1000); s=zeros(40, 1);
>> for n=1:40; tic(); [sumA A]= msum_oct(A); s(n) = toc(); end;
>> sum(s), mean(s), std(s)
>> for n=1:40; tic(); sumA = msum_oct(A); s(n) = toc(); end;
>> sum(s), mean(s), std(s)
</code></pre>
This shows a speed-up of a factor of approximately 2 on the computers I have tried. The actual speed-up depends on the array size, hardware, optimization flags, etc.   
</p>
</body>
</html>
