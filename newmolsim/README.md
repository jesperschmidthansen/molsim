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
