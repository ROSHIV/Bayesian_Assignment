%Task B

clc;
clear all;
close all;

% Task B _1,2
theta = [0.25,0.5];
gama = [0,1]
for i= 1:length(theta)
   b_g0 = Bernoulli_f(theta(i),gama(1)); 
   b_g1 = Bernoulli_f(theta(i),gama(2));
   b = [b_g0 b_g1];
   figure()
   stem(gama,b,'filled','red','linestyle','-.');
   ylabel('probebilities');
   xlabel('\gamma');
   title(strrep(['\theta=' num2str(theta(i)) ],' ','_'));
end 



% Task B _3
theta = 0:0.1:1; 
gama = [0,1];
for i= 1:length(theta)
   b_g0(i) = Bernoulli_f(theta(i),gama(1)); 
   b_g1(i) = Bernoulli_f(theta(i),gama(2));
end 
figure()
stem(theta,b_g0,'filled','red','linestyle','-.');
ylabel('liklihood');
xlabel('\theta');
hold on;
stem(theta,b_g1,'filled','green','linestyle','-.');
ylabel('liklihood'); 
xlabel('\theta');
legend('\gamma(0)','\gamma(1)');
hold off;

% Task B _4a,b,c
N = [10 1000,10000]
theta = 0.5;
for i =1:length(N)
d = rand(1,N(i));gama=zeros(1,N(i));
for j = 1:N(i)
    if d(j) > theta
       gama(j) = 1;
    end
end

for k= 1:length(gama)
   b(k) = Bernoulli_f(theta,gama(i)); 
   b_log(k)= log(theta).^gama(i).*log(1-theta).^(1-gama(i));
%    b_exp(i)= exp(theta).^gama(i).*exp(1-theta).^(1-gama(i));
end 
likelihood(i) = prod(b); 
log_likelihood(i) = sum(b_log);
% exp_likelihood = sum(b_exp);
fprintf('%d\n\n',N(i),likelihood(i),log_likelihood(i));
end

% Task B _4d
theta = 0:0.01:1;  
gama = {[1],[1,1],[1,1,0,1]}; % samples
figure()
for k = 1:size(gama,2)
   for i= 1:length(theta)
       z = sum(gama{k}); 
       N = length(gama{k}); % number of samples
       likelihood(i) = theta(i).^z.*(1-theta(i)).^(N-z);
   end
   
plot(theta,likelihood);
ylabel('likelihood');
xlabel('\theta');
hold on;
legend('\gamma[1]','\gamma[1,1]','\gamma[1,1,0,1]');
end



% Task B _5
a = 1;
b = 1;
theta_s = 0:0.001:1;
gama = {[1],[1,1],[1,1,0,1]}; % samples
figure()
for k = 1:size(gama,2)
    for i=1:length(theta_s)
    z = sum(gama{k});
    N = length(gama{k}); % number of samples
    log_posterior(i) = (z+a-1).*log(theta_s(i)) + (N-z+b-1).*log(1-theta_s(i)) - log(beta(z+a,N-z+b));
    end
    plot(theta_s,log_posterior);
    plot(theta_s,log_posterior);
    ylabel('log_posterior');
    xlabel('\theta_s');
    hold on;
    legend('\gamma[1]','\gamma[1,1]','\gamma[1,1,0,1]');
end

a=[250 18.25 1];
b=[250 6.75 1];
gama =[ones(17,1);zeros(3,1)];
N=length(gama);
z=sum(gama);
theta_s=0:0.1:1;
log_posterior = zeros(length(theta_s));
figure()
for j =1:length(a)
    
    for i = 1:length(theta_s)
    log_posterior(i) = (z+a(j)-1).*log(theta_s(i)) + (N-z+b(j)-1).*log(1-theta_s(i)) - log(beta(z+a(j),N-z+b(j)));
    end
    plot(theta_s,log_posterior);
    ylabel('log_posterior');
    xlabel('\theta_s');
    hold on;
    legend('N=20, z=17, a=b=250','N=20, z=17, a=18.25, b=6.75',' N=20,z=17, a=b=1');
    
end 
 


   
   

 



   
   