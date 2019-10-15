%Taska A

clc;
clear all;
close all;


keyset = {':(',':)','+|:(','+|:)','-|:(','-|:)'};
valueset = [0.001 0.999 0.99 0.05 1-0.99 1-0.05];

[p] = disease_test(keyset,valueset);  

p('+') = p('+|:)')*p(':)')+p('+|:(')*p(':('); %prob of positive test
p('-') = p('-|:)')*p(':)')+p('-|:(')*p(':(');

p(':)|+') = p('+|:)')*p(':)')/p('+'); % prob of having desease in terms of positive test 
p(':)|-') = p('-|:)')*p(':)')/p('-');

p(':(|+') = p('+|:(')*p(':(')/p('+');
p(':(|-') = p('-|:(')*p(':(')/p('-');


prior = [p(':('), p(':)')];
Test_number = 4;
T = round(rand(1,Test_number)); % T=0 (test negative), T=1 (test positive)
for i=1:Test_number
   
    if (T(i)==0) 
        p(':(|-')= p(':(|-')*p('-|:(')/(p('-|:(')*p(':(|-')+p('-|:)')*p(':)|-'));
        p(':)|-')= p(':)|-')*p('-|:)')/(p('-|:(')*p(':(|-')+p('-|:)')*p(':)|-'));
        healthy(i)=p(':(|-'); desease(i)=p(':)|-');
        
    else 
        p(':(|+')=p('+|:(')*p(':(|+')/(p('+|:(')*p(':(|+')+p('+|:)')*p(':)|+'));
        p(':)|+')=p('+|:)')*p(':)|+')/(p('+|:(')*p(':(|+')+p('+|:)')*p(':)|+'));
        healthy(i)=p(':(|+'); desease(i)=p(':)|+');
     
    end
end

 fprintf('%d\n\n\',T,healthy,desease);

