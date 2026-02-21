 <html>
<body>
<h1> molsim - Molecular dynamics with GNU Octave </h1>
<p>
<figure>
  <img src="resources/logo_0.png" style="width:40%" class="center">
</figure>
</p>

<h2>
 molsim supports simulations of
</h2>
<ul>
    <li>simple Lennard-Jones systems,</li>
    <li>molecular systems with bond, angle, and torsion potentials,</li>
    <li>confined flow systems, eg., Couette and Poiseuille flows,</li>
    <li>charged systems using shifted force,</li>
    <li>and more ...</li>
</ul>
</p>

<h2>Installation </h2>
<p>At the Octave prompt simply use the command </p>
<pre>
 <code>
  >> pkg install "https://github.com/jesperschmidthansen/molsim/archive/refs/tags/v&lt;version&gt;.tar.gz"
 </code> 
</pre>
<p>where &lt;version&gt; is the version number. 

<h2>Examples</h2>
Checkout the project example folder

<h2>Contribution</h2>
<p>
I encourage anyone who uses or plans to use molsim to submit problematic issues - this includes issues regarding the documentation. I also welcome contributions to the code for the project, whether it is core features or post simulation data analysis programs. 
</p>

<h2>Why MEX?</h2>
GNU Octave offers fantastic C++ interface with dynamically linked functions (DLDs). However, I find
the pure C MEX interface to produce faster running binaries. This is perhaps due to call-by-value and
call-by-reference strategies giving DLDs an additional copying overhead.    

Consider the two functions below. 
<table>
 <tr> <td> DLD</td><td>MEX</td></tr>
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

<p> dsds </p>


<h2>To-do</h2>
Octave now supports object oriented programming. molsim is under complete reconstructed to benefit
from this. Matlab compatibility is not a priority.

- [ ] Feature: Barostate
- [ ] Feature: Standard run time sample classes
- [ ] Feature: Electrostatic interactions using the Wolf scheme
- [ ] Feature: A set of molecular and atomic configurations 
- [ ] Feature: Molecular class for infrastructure
- [ ] Feature: DPD support (initiated)
- [ ] Revision: Class properties access. Should these be different from public?
- [ ] Revision: Define class constants with correct properties (Constant=true)
- [ ] Revision: All classes should have a disp method
- [ ] Revision: Consider whether methods should have specified properties
- [ ] Revision: Naming convensions (at the moment none)  
- [ ] Revision: ms_molconfig is a mess... 

</body>
</html>
