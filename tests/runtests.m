
dotests = [false,  
	   false,
	   false,
	   true];

itest = 1;
tests=["_test01.m"; "_test02.m"; "_test03.m"; "_test04.m"];

save test.mat itest dotests tests;

do 
   
   clear all;

   load test.mat
   
   if ( dotests(itest) )
     run(tests(itest,:));
   endif
   
   itest = 1 + itest;
   
   save test.mat itest dotests tests;
   
until ( itest==length(dotests)+1 )
  
