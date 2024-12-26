 <html>
<body>
<h1> molsim - Molecular dynamics with GNU Octave </h1>
<p>
<figure>
  <img src="doc/logo.png" alt="Trulli" style="width:30%">
</figure> 
</p>

<p>
molsim supports simulations of
</p>

<ul>
<li>simple Lennard-Jones systems,</li>
<li>molecular systems with bond, angle, and torsion potentials,</li>
<li>confined flow systems, eg., Couette and Poiseuille flows,</li>
<li>charged systems using shifted force and Wolf methods,</li>
<li>dissipative particle dynamics systems,</li>
<li>different ensembles,</li>
<li> and more . .</li>
</ul>

<h2>Installation </h2>
<p>At the Octave prompt simply use the command </p>
<pre>
 <code>
  >> pkg install "https://github.com/jesperschmidthansen/molsim/archive/refs/tags/v&lt;version&gt;.tar.gz"
 </code> 
</pre>
<p>where &lt;version&gt; is the version number. NOTE: Depdending on your system you may recieve warnings like
<pre>
 <code>
  note: expected 'const mwSize *' {aka 'const long long int *'} but argument is of type 'const long int *'
 </code> 
</pre>
These warnings are not always harmful and not a molsim issue. Proceed with care.

<h2>An example</h2>
An example of an NVE water simulation script

```octave
nloops = 1000; temp0 = 298.15/78.2;
cutoff= 2.5; sigma=1.0; epsilon=1.0; aw=1.0; cutoff_sf = 2.9;
lbond = 0.316; kspring = 68421; 
angle = 1.97; kangle = 490;

molsim('set', 'cutoff', cutoff_sf);
molsim('set', 'timestep', 0.0005);
molsim('set', 'exclusion', 'molecule'); 

molsim('set', 'omp', 4);

molsim('load', 'xyz', 'sys_water.xyz');  molsim('load', 'top', 'sys_water.top');

for n=1:nloops 
  molsim('reset')
  
  molsim('calcforce', 'lj', 'OO', cutoff, sigma, epsilon, aw);
  molsim('calcforce', 'coulomb', 'sf', cutoff_sf);
  molsim('calcforce', 'bond', 0, lbond, kspring);
  molsim('calcforce', 'angle', 0, angle, kangle);
  
  molsim('integrate', 'leapfrog');
end

molsim('clear');
```
<p> <b>IMPORTANT NOTE</b>: The 'sys_water.xyz' configuration file and 'sys_water.top' topology file must be in
same directory from where you execute the script. They can be found under the project's resource/ folder </p>
<p> For further explanation check out the package tutorial under the project's doc/ folder </p> 

<h2>Contribution</h2>
<p>
I encourage anyone who uses or plans to use molsim to submit problematic issues - this includes issues regarding the documentation. I also welcome contributions to the code for the project, whether it is core features (seplib), post simulation data analysis programs, or extending the molsim wrapper. 
</p>

<h2>To-do</h2>
Octave now supports object oriented programming. molsim is under complete reconstructed to benefit from this, see folder newmolsim. Matlab compability
is relaxed.

- [ ] Feature: Molecular force fields
- [ ] Feature: Barostate
- [ ] Feature: Standard run time sample classes
- [ ] Feature: Electrostatics
- [X] Revision: endX -> end
- [ ] Revision: Class properties access. Should these be different from public?
- [ ] Test/examples: Example folder with different scripts





</body>
</html>
