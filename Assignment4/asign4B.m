% Asignement4
% Task B

clc;
close all;
clear all;

%input parameter
gama_y = [1,0,1,1,0,1,1,1,0,1,1,1,1,1];% data 1 
gama_z = [1,0,0,0,0,0,0,1,1,0]; % data 2
a = 1; b =1;
N_samples=5000;
z_y = sum(gama_y); 
z_z = sum(gama_z);   
N_data_y = length(gama_y); % number of samples
N_data_z = length(gama_z);

% function
Beta_posterior_f_y = @(theta,a,b,z,N_data) (theta.^(z+a-1).*((1-theta).^(N_data-z+b-1)))./(beta(z+a,N_data-z+b));
Beta_posterior_f_z = @(theta,a,b,z_z,N_data_z) (theta.^(z_z+a-1).*((1-theta).^(N_data_z-z_z+b-1)))./(beta(z_z+a,N_data_z-z_z+b));

theta_ini=0.5;
theta_y = slicesample(theta_ini,N_samples,'burnin', 1000,'pdf',@(theta) Beta_posterior_f_y(theta,a,b,z_y,N_data_y),'width',0.01);
theta_z = slicesample(theta_ini,N_samples,'burnin', 1000,'pdf',@(theta) Beta_posterior_f_z(theta,a,b,z_z,N_data_z),'width',0.01);
nbins=100;
h=histogram(theta_y,nbins,'normalization','pdf');
hold on
h=histogram(theta_z,nbins,'normalization','pdf');


% probrbility of theta_y>0.5
p_y = length(find(theta_y>=0.5))/length(theta_y);
p_z = length(find(theta_z>=0.5))/length(theta_z);

% HDI ( highst density interval)
SEM = std(theta_y)/sqrt(length(theta_y));        % Standard Error
ts = tinv([0.005  0.95],length(theta_y)-1);      % T-Score
CI_y = mean(theta_y) + ts*SEM; 
plot([0,1],CI_y,'color','r','LineWidth',2,'LineStyle',':','Color','b');
SEM = std(theta_z)/sqrt(length(theta_z));               % Standard Error
ts = tinv([0.005  0.95],length(theta_z)-1);      % T-Score
CI_z = mean(theta_z) + ts*SEM; 
plot([0,1],CI_z,'color','r','LineWidth',2,'LineStyle',':','Color','r');
xlabel('\theta');
grid on;
legend('\theta_y','\theta_z');

% d_theta / (theta_y)-(theta_z)
figure()
d_theta = theta_y-theta_z;
h=histogram(d_theta,nbins,'normalization','pdf');
xlabel('\theta');


% HDI ( highst density interval)
SEM = std(d_theta)/sqrt(length(d_theta));        
ts = tinv([0.005  0.95],length(d_theta)-1);      
CI_d = mean(d_theta) + ts*SEM; 
hold on 
plot([0,1],CI_d,'color','r','LineWidth',2,'LineStyle',':','Color','b');
  
  


