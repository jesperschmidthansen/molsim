#
# pre_install function 
#

function pre_install(dummy)

	run("src/mhtxt.m");
	system("mv src/*.m inst/");

end

