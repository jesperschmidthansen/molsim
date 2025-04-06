#
# pre_install function 
#

function pre_install(dummy)

	run("src/mhtxt.m");
	system("mv src/ms_*.m inst/");
	system("mv src/mhtxt.m inst/");	
end

