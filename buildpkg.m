
tar("molsim.tar", "molsim");
gzip("molsim.tar");
pkg build ./ molsim.tar.gz
gunzip("molsim.tar.gz");
untar("molsim.tar");
