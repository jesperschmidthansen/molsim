
dotests = [true,true,true,true,true,true];

itest = 1;
tests=["_test01.m"; "_test02.m"; "_test03.m"; "_test04.m" ; "_test05.m"; "_test06.m"];

save test.mat itest dotests tests;

do 
   
   clear all;

   load test.mat
   
   if ( dotests(itest) )
     run(tests(itest,:));
   endif

   load test.mat
    
   itest = 1 + itest;
   
   save test.mat itest dotests tests;
   
until ( itest==length(dotests)+1 )
  
