% Asignement4
% Task A

clc;
close all;
clear all;

%input parameter
gama =[ones(17,1);zeros(3,1)]; %data
N_data = length(gama); % number of data
z = sum (gama); 
a = 100; b =100;
N_samples=5000; std=0.1;

%plot prior(beta),posterior(beta),likelihood 
theta=0:0.001:1;

%function 
Beta_f = @(theta,a,b) ((theta.^(a-1)).*((1-theta).^(b-1)))./beta(a,b);
Bernoulli_liklihood_f = @(theta,N_data,z) (theta).^z.*((1-theta)).^(N_data-z);
Beta_posterior_f = @(theta,a,b,z,N_data) (theta.^(z+a-1).*((1-theta).^(N_data-z+b-1)))./(beta(z+a,N_data-z+b));


figure()
for j = 1:3
    
    for i = 1:length(theta)
        prior(i) = Beta_f(theta(i),a,b);
        likelihood(i) = Bernoulli_liklihood_f(theta(i),N_data,z);
        posterior(i) = Beta_posterior_f(theta(i),a,b,z,N_data);
    end
    
end

subplot(3,1,1);
plot(theta,prior);
ylabel(strrep(['p(D|\theta' num2str(a) ',' num2str(b) ')' ],' ','_'));
xlabel('\theta');
title ('Prior(beta)');
legend(strrep(['N=' num2str(N_data) ',' 'z=' num2str(z) ')' ],' ','_'));
hold on
subplot(3,1,2);
plot(theta,likelihood);
ylabel('p(D|\theta)');
xlabel('\theta');
title ('Likelihood(Bernoulli)');
hold on 
subplot(3,1,3);
plot(theta,posterior,'Color', 'r');
ylabel(strrep(['dbeta(\theta|' num2str(z+a) ',' num2str(N_data-z+b) ')' ],' ','_'));
xlabel('\theta');
title ('Posterior(beta)');
hold on 

%metropolis_sampling

% theta_ini=0.5;
% temp=1;
% 
% pdf_prior = makedist('Beta','a',a,'b',b);
% posterior_sample(1)=theta_ini;
% 
% for i=1:N_samples 
%     d_theta=std*randn(1);
%     theta_p=theta_ini+d_theta;
%     if theta_p>=0&&theta_p<=1
%        alpha=(pdf(pdf_prior,theta_p)*Bernoulli_liklihood_f(theta_p,N_data,z))/(pdf(pdf_prior,theta_ini)*Bernoulli_liklihood_f(theta_ini,N_data,z));
%        Jump(temp)=min(1,alpha);
%        sample=rand(1);
%        if sample <= Jump(temp)
%        temp=temp+1;
%        posterior_sample(temp)=theta_p;
%        theta_ini=theta_p;
%         end
%     end
% end
% 
% nbins=100;
% h=histogram(posterior_sample,nbins,'normalization','pdf');
% h.FaceColor = 'b';
% h.EdgeColor = 'g'; 
% legend('posterior(beta)','posterior_sample');
% theta_ini=0.5;N_samples=1000;


% slice smapling (Matlab function)
theta_ini=0.5;
theta = slicesample(theta_ini,N_samples,'burnin', 1000,'pdf',@(theta) Beta_posterior_f(theta,a,b,z,N_data),'width',0.01);
nbins=100;
h=histogram(theta,nbins,'normalization','pdf');
h.FaceColor = 'b';
h.EdgeColor = 'g'; 
legend('posterior(beta)','posterior-sample');



