% Asignement6
clc;
close all;
clear all;

%input parameter
Y=[607	583	521	494	369	782	570	678	467	620	425	395	346	361	310	300	382	294	315	323	421	339	398	328	335	291	329	310	294	321	286	349	279	268	293	310	259	241	243	272	247	275	220	245	268	357	273	301	322	276	401	368	149	507	411	362	358	355	362	324	332	268	259	274	248	254	242	286	276	237	259	251	239	247	260	237	206	242	361	267	245	331	357	284	263	244	317	225	254	253	251	314	239	248	250	200	256	233	427	391	331	395	337	392	352	381	330	368	381	316	335	316	302	375	361	330	351	186	221	278	244	218	126	269	238	194	384	154	555	387	317	365	357	390	320	316	297	354	266	279	327	285	258	267	226	237	264	510	490	458	425	522	927	555	550	516	548	560	545	633	496	498	223	222	309	244	207	258	255	281	258	226	257	263	266	238	249	340	247	216	241	239	226	273	235	251	290	473	416	451	475	406	349	401	334	446	401	252	266	210	228	250	265	236	289	244	327	274	223	327	307	338	345	381	369	445	296	303	326	321	309	307	319	288	299	284	278	310	282	275	372	295	306	303	285	316	294	284	324	264	278	369	254	306	237	439	287	285	261	299	311	265	292	282	271	268	270	259	269	249	261	425	291	291	441	222	347	244	232	272	264	190	219	317	232	256	185	210	213	202	226	250	238	252	233	221	220	287	267	264	273	304	294	236	200	219	276	287	365	438	420	396	359	405	397	383	360	387	429	358	459	371	368	452	358	371];
ind = [1	1	1	1	1	2	2	2	2	2	3	3	3	3	3	3	3	3	3	4	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	6	6	6	6	6	6	6	6	6	7	7	7	7	7	8	8	8	8	8	9	9	9	9	9	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	13	13	13	13	13	13	13	13	13	13	14	14	14	14	14	14	14	14	14	14	14	14	14	15	15	15	15	15	15	16	16	16	16	16	17	17	17	17	17	18	18	18	18	18	19	19	19	19	19	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	21	21	21	21	21	21	21	21	21	21	22	22	22	22	22	22	22	22	22	22	22	22	23	23	23	23	23	23	23	23	23	23	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	25	25	25	25	25	25	25	25	25	25	25	25	25	26	26	26	26	26	27	27	27	27	27	28	28	28	28	28	28	28	28	28	29	29	29	29	29	29	29	29	29	29	30	30	30	30	30	30	31	31	31	31	31	32	32	32	32	32	33	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34]; 
child_j = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
child_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


%% Function 

%Prior
% log_p_Sigma =  @(Sigma) log(real(Sigma>0));
% log_p_tau = @(tau) 0;
log_p_mu = @(mu) 0;
log_p_phi = @(phi) log(real(phi>0)); %ok<ELARLOG>
log_p_theta = @(theta,mu,tau) sum(-log(tau)-0.5.*((theta-mu)./tau).^2);
log_p_theta_classified = @(theta,mu,tau,phi) sum(-log(tau)-0.5.*((theta-(mu+phi.*child_j))./tau).^2);
%log posterior (hirarchical model)
log_posterior = @(bita) log_p_likelihood(bita(1),bita(4:end)) + log_p_theta(bita(4:end),bita(3),bita(2)) + log_p_mu(bita(3)) + log_p_Sigma(bita(1)) +log_p_tau(bita(2));
%log posterior (classified)
log_posterior_classified = @(bita) log_p_likelihood(bita(1),bita(5:end)) + log_p_theta_classified(bita(5:end),bita(3),bita(2),bita(4)) + log_p_mu(bita(3)) + log_p_Sigma(bita(1)) +log_p_tau(bita(2)) + log_p_phi(bita(4));


%% slice sampling_classified  function slice_sample.m (Adult and Kids)

%Parameters 
sigma_ini_c = 0.2;
tau_ini_c = 0.2;
mu_ini_c = 6;
phi_ini_c = 0.6;
theta_ini_c = 6;
bita_ini_classified = [sigma_ini_c,tau_ini_c,mu_ini_c,phi_ini_c,theta_ini_c.*ones(1,34)];
N_samples = 10000;
widths = 1; burnIn=0;

rnd_classified = slice_sample(N_samples,burnIn,log_posterior_classified,bita_ini_classified,widths,false); 

sigma_s_c = rnd_classified(1,:);
tau_s_c = rnd_classified(2,:);
mu_s_c = rnd_classified(3,:);
phi_s_c = rnd_classified(4,:);
theta_s_c = rnd_classified(5:end,:);



%% slice sampling  function slice_sample.m 

%Parameters 
sigma_ini = 0.2;
tau_ini = 0.2;
mu_ini = 6;
theta_ini = 6;
bita_ini = [sigma_ini,tau_ini,mu_ini,theta_ini.*ones(1,34)];
N_samples = 10000;
widths = 1; burnIn=0;

rnd = slice_sample(N_samples,burnIn,log_posterior,bita_ini,widths,false); 

sigma_s = rnd(1,:);
tau_s = rnd(2,:);
mu_s = rnd(3,:);
theta_s = rnd(4:end,:);




%% Task1

%plot
Mean_phi_s_c = mean(phi_s_c);
Median_phi_s_c = median(phi_s_c);
figure()
nbins=500;
h=histogram(phi_s_c,nbins,'Normalization','pdf');
xlabel('\phi','FontSize', 18);
hold on
xline((Mean_phi_s_c),'r',strrep(['Mean =' num2str(Mean_phi_s_c) ],' ','_'), 'LineWidth', 2, 'LabelOrientation', 'horizontal');
hold on
xline(Median_phi_s_c,'-- b',strrep(['Median =' num2str(Mean_phi_s_c) ],' ','_'), 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
legend('\phi');
hold on
HDI_1 = prctile(phi_s_c,97.5);
HDI_2 = prctile(phi_s_c,2.5);xline(HDI_1,'LineWidth', 2);xline(HDI_2, 'LineWidth', 2);
hold on 


%% Task2

Mean_tau_s = mean(tau_s);
Median_tau_s = median(tau_s);
figure()
nbins=500;
h=histogram(tau_s,nbins,'Normalization','pdf');
xlabel('\tau','FontSize', 14);
hold on
xline((Mean_tau_s),'r',strrep(['Mean =' num2str(Mean_tau_s) ],' ','_'), 'LineWidth', 2, 'LabelOrientation', 'horizontal');
hold on
xline(Median_tau_s,'-- b',strrep(['Median =' num2str(Median_tau_s) ],' ','_'), 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
legend('\tau-asi5','Mean','Median');
hold on
HDI_1 = prctile(tau_s,97.5);
HDI_2 = prctile(tau_s,2.5);xline(HDI_1,'');xline(HDI_2,'');
hold on

Mean_tau_s_c = mean(tau_s_c);
Median_tau_s_c = median(tau_s_c);
figure()
nbins=500;
h=histogram(tau_s_c,nbins,'Normalization','pdf');
xlabel('\tau','FontSize', 14);
hold on
xline((Mean_tau_s_c),'r',strrep(['Mean =' num2str(Mean_tau_s_c) ],' ','_'), 'LineWidth', 2, 'LabelOrientation', 'horizontal');
hold on
xline(Median_tau_s_c,'-- b',strrep(['Median =' num2str(Mean_phi_s_c) ],' ','_'), 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
legend('\tau-asi6','Mean','Median');
hold on
HDI_1 = prctile(tau_s_c,97.5);
HDI_2 = prctile(tau_s_c,2.5);xline(HDI_1,'');xline(HDI_2,'');
hold on


%% Task3
theta_n = mu_s+randn(1,N_samples).*tau_s; % done  in assignmetn 5
theta_adult = mu_s_c+randn(1,N_samples).*tau_s_c;
theta_kids = mu_s_c+phi_s_c +randn(1,N_samples).*tau_s_c;

Mean_theta_n = mean(theta_n);
Median_theta_n = median(theta_n);
figure()
nbins=100;
h=histogram(theta_n,nbins,'Normalization','pdf');
xlabel('\theta','FontSize', 14);
xline((Mean_theta_n),'r', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(Median_theta_n,'-- b', 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
legend('\theta-asi5', 'Mean', 'Median');
HDI_1 = prctile(theta_n,97.5);
HDI_2 = prctile(theta_n,2.5);
%xline(HDI_1,'');xline(HDI_2,'');
hold on 

Mean_theta_adult = mean(theta_adult);
Median_theta_adult = median(theta_adult);
h=histogram(theta_adult,nbins,'Normalization','pdf');
xlabel('\theta','FontSize', 14);
xline((Mean_theta_adult),'r', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(Median_theta_adult,'-- b', 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
HDI_1 = prctile(theta_adult,97.5);
HDI_2 = prctile(theta_adult,2.5);
%xline(HDI_1,'');xline(HDI_2,'');
hold on 

Mean_theta_kids = mean(theta_kids);
Median_theta_kids = median(theta_kids);
h=histogram(theta_kids,nbins,'Normalization','pdf');
xlabel('\theta','FontSize', 14);
xline((Mean_theta_kids),'--r', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(Median_theta_kids,'-- b', 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
HDI_1 = prctile(theta_kids,97.5);
HDI_2 = prctile(theta_kids,2.5);
%xline(HDI_1,'');xline(HDI_2,'');
hold on 

%% Task4
y_adult = exp(mu_s_c+randn(1,N_samples).*tau_s_c+randn(1,N_samples).*sigma_s_c);
y_kids = exp(mu_s_c+phi_s_c+randn(1,N_samples).*tau_s_c+randn(1,N_samples).*sigma_s_c);
Mean_y_kids = mean(y_kids);
Median_y_kids = median(y_kids);
figure()
nbins=500;
h=histogram(y_kids,nbins,'Normalization','pdf');
xlabel('\tau','FontSize', 14);
hold on
xline((Mean_y_kids),'r', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(Median_y_kids,'-- b', 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
%legend('\adult', 'Mean', 'Median');
HDI_1 = prctile(theta_n,97.5);
HDI_2 = prctile(theta_n,2.5);
%xline(HDI_1,'');xline(HDI_2,'');
legend('\theta-adult');
hold on

Mean_y_adult = mean(y_adult);
Median_y_adult = median(y_adult);
h=histogram(y_adult,nbins,'Normalization','pdf');
xlabel('\tau','FontSize', 14);
hold on
xline((Mean_y_adult),'r', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
hold on
xline(Median_y_adult,'-- b', 'LineWidth', 2,'LabelOrientation', 'horizontal','LabelHorizontalAlignment', 'left' );
HDI_1 = prctile(y_adult,97.5);
HDI_2 = prctile(y_adult,2.5);
%xline(HDI_1,'');xline(HDI_2,'');
hold on 



%% log_liklihod
function [p] = log_p_likelihood(Sigma,theta)

%input parameter
Y=[607	583	521	494	369	782	570	678	467	620	425	395	346	361	310	300	382	294	315	323	421	339	398	328	335	291	329	310	294	321	286	349	279	268	293	310	259	241	243	272	247	275	220	245	268	357	273	301	322	276	401	368	149	507	411	362	358	355	362	324	332	268	259	274	248	254	242	286	276	237	259	251	239	247	260	237	206	242	361	267	245	331	357	284	263	244	317	225	254	253	251	314	239	248	250	200	256	233	427	391	331	395	337	392	352	381	330	368	381	316	335	316	302	375	361	330	351	186	221	278	244	218	126	269	238	194	384	154	555	387	317	365	357	390	320	316	297	354	266	279	327	285	258	267	226	237	264	510	490	458	425	522	927	555	550	516	548	560	545	633	496	498	223	222	309	244	207	258	255	281	258	226	257	263	266	238	249	340	247	216	241	239	226	273	235	251	290	473	416	451	475	406	349	401	334	446	401	252	266	210	228	250	265	236	289	244	327	274	223	327	307	338	345	381	369	445	296	303	326	321	309	307	319	288	299	284	278	310	282	275	372	295	306	303	285	316	294	284	324	264	278	369	254	306	237	439	287	285	261	299	311	265	292	282	271	268	270	259	269	249	261	425	291	291	441	222	347	244	232	272	264	190	219	317	232	256	185	210	213	202	226	250	238	252	233	221	220	287	267	264	273	304	294	236	200	219	276	287	365	438	420	396	359	405	397	383	360	387	429	358	459	371	368	452	358	371];
ind = [1	1	1	1	1	2	2	2	2	2	3	3	3	3	3	3	3	3	3	4	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	6	6	6	6	6	6	6	6	6	7	7	7	7	7	8	8	8	8	8	9	9	9	9	9	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	13	13	13	13	13	13	13	13	13	13	14	14	14	14	14	14	14	14	14	14	14	14	14	15	15	15	15	15	15	16	16	16	16	16	17	17	17	17	17	18	18	18	18	18	19	19	19	19	19	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	21	21	21	21	21	21	21	21	21	21	22	22	22	22	22	22	22	22	22	22	22	22	23	23	23	23	23	23	23	23	23	23	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	25	25	25	25	25	25	25	25	25	25	25	25	25	26	26	26	26	26	27	27	27	27	27	28	28	28	28	28	28	28	28	28	29	29	29	29	29	29	29	29	29	29	30	30	30	30	30	30	31	31	31	31	31	32	32	32	32	32	33	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34];
% ZLog_Y=  (log(Y)-mean(log(Y)))./std(log(Y)); % standardize log_Y
ZLog_Y=  log(Y); % log_Y
p = 0;
for i=1:length(Y)
    j=ind(i);
    p = p+(-log(Sigma)-0.5.*((ZLog_Y(i)-theta(j))./Sigma)^2); % Nomal probebility density function (liklihood)
end
end

%% log_sigma
function output = log_p_Sigma(Sigma)
if Sigma>0
    output = 0;
else
    output = -inf;
end
end

%% log_tau
function output = log_p_tau(tau)
if tau>0
    output = 0;
else
    output = -inf;
end
end

