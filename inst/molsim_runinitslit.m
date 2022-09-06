
clear all;


function [npart str] = readheader(fin, opt)

  npart = fscanf(fin, "%d\n", "C");

  if strcmp(opt, "string")
    str = fgets(fin, 1024);
  else 
    str = '\0';
  endif

endfunction

function add_wallforce(lbox)
  
  x = molsim('get', 'positions');
  
  f = 1./x(:,3).^12 - 1./(x(:,3)-lbox).^12;
  
  molsim('add', 'force', f, 3);
  
endfunction


function write_config(lbox, boff)
  
  str=sprintf("sep_lattice -n=%d,%d,3 -l=%f,%f,3.0 -f=wall.xyz", ...
		int32(lbox(1)), int32(lbox(2)), lbox(1), lbox(2));

  system(str);

  fin_mol = fopen("molecules.xyz", 'r');
  fin_wall = fopen("wall.xyz", 'r');

  npart_wall = readheader(fin_wall, 'string');
  npart_mol = readheader(fin_mol, 'string');

  fout = fopen("slitpore.xyz", 'w');

  fprintf(fout, "%d\n", npart_wall*2+npart_mol);
  fprintf(fout, "%f %f %f\n", lbox(1), lbox(2), lbox(3)*2);

  offset = 2.0 + boff; zmax = 0.0;
  for n=1:npart_mol
    [t, x, y, z, vx, vy, vz, m, q] = fscanf(fin_mol, "%c %f %f %f %f %f %f %f %f\n", "C");
    fprintf(fout, "%c %f %f %f %f %f %f %f %f\n",t, x, y, z+offset, vx, vy, vz, m, q);

    if z > zmax
      zmax = z;
    end
    
  endfor

  offset = zmax + 2*boff + 2.0;

  for n=1:npart_wall
    [t, x, y, z, vx, vy, vz, m, q] = fscanf(fin_wall, "%c %f %f %f %f %f %f %f %f\n", "C");
    fprintf(fout, "w %f %f %f %f %f %f %f %f\n", x, y, z, vx, vy, vz, m, q);
    fprintf(fout, "W %f %f %f %f %f %f %f %f\n", x, y, offset+z, vx, vy, vz, m, q)
  endfor
  
  fclose(fin_mol);fclose(fin_wall);fclose(fout);

endfunction


function Lbox = compress(dens0, lbox_z, xyzfile, topfile)

  temp0 = 4.0;

  lb = 0.4; ks = 2000.0;
  angle = 1.9; kangle = 400.0;
  torsionparam = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070];

  molsim('set', 'molconfig', xyzfile, topfile, 1000, 0.01, 42);
  molsim('load', 'xyz', 'start.xyz');
  molsim('load', 'top', 'start.top');

  molsim('set','timestep', 0.001);
  molsim('set', 'temperature', temp0);
  molsim('set', 'exclusion', 'molecule');
  molsim('set', 'compressionfactor', 0.99995);

  npart = molsim('get', 'numbpart');

  Lbox = molsim('get', 'box');
  lbox_xy = sqrt(npart/(lbox_z*dens0)); 
  
  while ( Lbox(1) > lbox_xy || Lbox(3) > lbox_z )

    molsim('reset')
    
    molsim('calcforce', 'lj', 'CC', 2.5, 1.0, 1.0, 1.0);
    
    molsim('calcforce', 'bond', 0, lb, ks);
    molsim('calcforce', 'angle', 0, angle, kangle);
    molsim('calcforce', 'torsion', 0, torsionparam);
    
    add_wallforce(molsim('get', 'box')(3));
    
    molsim('integrate', 'leapfrog')
    molsim('thermostat', 'relax', 'C', temp0, 0.001);
    
    molsim('compress', lbox_xy, 1);
    molsim('compress', lbox_xy, 2);
    molsim('compress', lbox_z, 3);

    Lbox = molsim('get', 'box');
  end

  molsim('save', 'C', 'molecules.xyz');
  molsim('clear');
  
end


densFluid = 1.18;
lengthWidth = 15.0;
xyzFile = 'butane.xyz';
topFile = 'butane.top';
bodyCentered = false;


lbox = compress(densFluid, lengthWidth, xyzFile, topFile);
write_config(lbox,1.2)

