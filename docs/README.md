 <html>
<body>

UNDER DEVELOPMENT 

This text is not meant to introduce molecular dynamics; such introductions can be found in many standard books. In brief, the basic idea is to solve the classical 
equation of motion of an ensemble of interacting particles. In the simplest form this means solving (numerically) Newton's second law

$$
   \frac{\mathrm{d}\mathbf{r}_i}{\mathrm{d} t} = \mathbf{v}_i \ , \ \
   \frac{\mathrm{d}\mathbf{p}_i}{\mathrm{d} t} = \mathbf{f}_i \ ,
$$

where $\mathbf{r}_i, \mathbf{v}_i, \mathbf{p}_i$ and $\mathbf{f}_i$ are the particle position, velocity, momentum and force acting on the particle, respectively. In a
standard simulation we solve this set of differential equations by (i) evaluating the forces acting on the particles, and (ii) from this integrate forward in time. The following pseudo code lists the basic idea

<pre><code>
Set simulation parameters
Set initial configuration position r and momenta p
 
do (as many times as we want)
   f <- calcforce(r)
   r, p <- integrate(f,p)
done
</pre></code>
The force is given by the gradient of the potential energy function $U$ by $\mathbf{f} = - \nabla U$. In molsim the energy function is 

$$
 U(\mathbf{r}_i, r_{ij}, \ldots) =  U_\mathrm{lattice} + U_\mathrm{vWaals} + U_{\mathrm{coulomb}} + U_\mathrm{bonds} + U_\mathrm{angles} +  U_\mathrm{torsion}
$$

The table shows the terms 
<table>
 <tr> <td> Potential function </td> <td> User supplied parameters </td> <td> molsim function </td> </tr>
 <tr>
  <td> $U_\mathrm{lattice} = \sum_\mathrm{sites} \frac{1}{2}k_0 (\mathbf{r}_i - \mathbf{r}_\text{lattice})^2$  </td>
 <td> $k_0$</td>
 <td> molsim.atoms.tether(atom type, $k_0$) </td>
 </tr>
 <tr>
  <td>$U_\mathrm{vWaals} =  \sum_{i,j \, \mathrm{pairs}} 4\epsilon\left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - a_w \left(\frac{\sigma}{r_{ij}}\right)^{6}\right]$ </td>
  <td> $r_\text{cutoff} $, $\epsilon$, $\sigma$, $a_w$</td>
  <td> molsim.lennardjones(atoms types, [$r_\text{cutoff}$, $\epsilon$, $\sigma$, $a_w$]); </td>
 </tr> 
<tr> 
 <td> $U_{\mathrm{coulomb}} = \sum_{i,j \, \mathrm{pairs}}\frac{q_iq_j}{r_{ij}}$</td> 
<td> $r_\text{cutoff}$ </td>
<td> molsim.sfcoulomb($r_\text{cutoff}$) </td>
</tr>
<tr><td> $U_{\mathrm{bonds}} =\sum_{\mathrm{bonds}} \frac{1}{2} k_{s}(r_{ij} - l_0)^2$ </td>
<td> $k_s$, $l_0$ </td>
<td> molsim.harmonicbond()</td>
</tr>
<tr> <td> $U_{\mathrm{angles}}=\frac{1}{2}\sum_{\mathrm{angles}} k_{\theta} (\cos(\theta) - \cos(\theta_0))^2$ </td>
<td> $k_\theta$, $\theta_0$ </td>
<td> molsim.harmonicangle() </td>
</tr>
 
<tr> 
 <td> $U_\mathrm{torsion}=\sum_{\mathrm{angles}} \sum_{n=0}^5 c_n \cos^n(\pi-\phi)$ </td>
<td> $c_n$ </td></tr>
<td>  sim.ryckbell();</td>
</table>



 
 






<h2>Why MEX?</h2>
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
