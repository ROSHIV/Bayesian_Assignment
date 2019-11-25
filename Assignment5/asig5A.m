% Asignement5
clc;
close all;
clear all;

%input parameter
Y=[607	583	521	494	369	782	570	678	467	620	425	395	346	361	310	300	382	294	315	323	421	339	398	328	335	291	329	310	294	321	286	349	279	268	293	310	259	241	243	272	247	275	220	245	268	357	273	301	322	276	401	368	149	507	411	362	358	355	362	324	332	268	259	274	248	254	242	286	276	237	259	251	239	247	260	237	206	242	361	267	245	331	357	284	263	244	317	225	254	253	251	314	239	248	250	200	256	233	427	391	331	395	337	392	352	381	330	368	381	316	335	316	302	375	361	330	351	186	221	278	244	218	126	269	238	194	384	154	555	387	317	365	357	390	320	316	297	354	266	279	327	285	258	267	226	237	264	510	490	458	425	522	927	555	550	516	548	560	545	633	496	498	223	222	309	244	207	258	255	281	258	226	257	263	266	238	249	340	247	216	241	239	226	273	235	251	290	473	416	451	475	406	349	401	334	446	401	252	266	210	228	250	265	236	289	244	327	274	223	327	307	338	345	381	369	445	296	303	326	321	309	307	319	288	299	284	278	310	282	275	372	295	306	303	285	316	294	284	324	264	278	369	254	306	237	439	287	285	261	299	311	265	292	282	271	268	270	259	269	249	261	425	291	291	441	222	347	244	232	272	264	190	219	317	232	256	185	210	213	202	226	250	238	252	233	221	220	287	267	264	273	304	294	236	200	219	276	287	365	438	420	396	359	405	397	383	360	387	429	358	459	371	368	452	358	371];
ind = [1	1	1	1	1	2	2	2	2	2	3	3	3	3	3	3	3	3	3	4	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	5	6	6	6	6	6	6	6	6	6	7	7	7	7	7	8	8	8	8	8	9	9	9	9	9	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	10	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	11	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	12	13	13	13	13	13	13	13	13	13	13	14	14	14	14	14	14	14	14	14	14	14	14	14	15	15	15	15	15	15	16	16	16	16	16	17	17	17	17	17	18	18	18	18	18	19	19	19	19	19	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	20	21	21	21	21	21	21	21	21	21	21	22	22	22	22	22	22	22	22	22	22	22	22	23	23	23	23	23	23	23	23	23	23	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	24	25	25	25	25	25	25	25	25	25	25	25	25	25	26	26	26	26	26	27	27	27	27	27	28	28	28	28	28	28	28	28	28	29	29	29	29	29	29	29	29	29	29	30	30	30	30	30	30	31	31	31	31	31	32	32	32	32	32	33	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34]; 


%% Function 

%Prior
log_p_Sigma =  @(Sigma) 0;
log_p_mu = @(mu) 0;
log_p_tau = @(tau) 0;
log_p_theta = @(theta,mu,tau) sum(-log(tau)-0.5.*((theta-mu)./tau).^2);
%log posterior
log_posterior = @(bita) log_p_likelihood(bita(1),bita(4:end)) + log_p_theta(bita(4:end),bita(3),bita(2)) + log_p_mu(bita(3)) + log_p_Sigma(bita(1)) +log_p_tau(bita(2));


%% slice sampling.. function slice_sample.m 

%Parameters 
sigma_ini = 0.2;
tau_ini = 0.2;
mu_ini = 6;
theta_ini = 6;
bita_ini = [sigma_ini,tau_ini,mu_ini,theta_ini.*ones(1,34)];
N_samples = 1500;
widths = 1; burnIn=0;

rnd = slice_sample(N_samples,burnIn,log_posterior,bita_ini,widths,false); 

sigma_s = rnd(1,:);
tau_s = rnd(2,:);
mu_s = rnd(3,:);
theta_s = rnd(4:end,:);




%% TaskA1

%Mode, Mean, HDI(highest Desity interval)
dude = exp(theta_s(4,:)+0.5.*sigma_s.^2);
Mean_dude= mean(dude);
% Mode_dude= mode(dude,'all');
Median_dude= median(dude);

%plot
figure()
nbins=100;
h=histogram(dude,nbins,'Normalization','pdf');
title('Expected reaction time for individual 4','FontSize', 14);
ylabel('Density','FontSize',14);
xlabel('Expected reacation time(ms), exp(\theta4-\sigma^2/2)','FontSize', 14);
hold on
xline(Mean_dude,'--r');
hold on
xline(Median_dude,'-- b');
legend('posterior','Mean','Median');
hold on
HDI_1=prctile(dude,97.5);
HDI_2=prctile(dude,2.5);xline(HDI_1,'');xline(HDI_2,'');


%% TaskA2a

% TaskA2 group reaction time 
group = exp(mu_s+(0.5.*(tau_s.^2))+0.5.*(sigma_s.^2));
Median_group = median(group);
Mean_group= mean(group);
figure()
nbins=100;
h=histogram(group,nbins,'Normalization','pdf');
title('Expected reaction time for group','FontSize', 14);
ylabel('Density','FontSize',14);
xlabel('Expected reacation time(ms), exp(\mu+\tau.^2/2+\sigma^2/2) ','FontSize', 14);
hold on
xline(Mean_group,'--r');
hold on
xline(Median_group,'-- b');
legend('posterior','Mean','Median');
hold on
HDI_1=prctile(group,97.5);
HDI_2=prctile(group,2.5);xline(HDI_1,'');xline(HDI_2,'');

% TaskA2ai Random Individual ( expected reaction timee)

theta_n = mu_s+randn(1,N_samples).*tau_s;
dude_random = exp(theta_n +(sigma_s.^2)/2); 
Median_dude_random = median(dude_random);
Mean_dude_random = mean(dude_random);
figure()
nbins=100;
h=histogram(dude_random,nbins,'Normalization','pdf');
title('Expected reaction time for random individual','FontSize', 14);
ylabel('Density','FontSize',14);
xlabel('Expected reacation time(ms)','FontSize', 14);
hold on
xline(Mean_dude_random,'--r');
hold on
xline(Median_dude_random,'-- b');
legend('posterior','Mean','Median');
hold on
HDI_1=prctile(dude_random,97.5);
HDI_2=prctile(dude_random,2.5);xline(HDI_1,'');xline(HDI_2,'');


% TaskA2aii Random Individual (predicted reaction timee)

% pred_random = (theta_n +randn(1,N_samples).*tau_s); 
% Median_dude_random = median(pred_random);
% Mean_dude_random = mean(pred_random);
% figure()
% nbins=100;
% h=histogram(pred_random,nbins,'Normalization','pdf');
% title('Expected reaction time for random individual','FontSize', 14);
% ylabel('Density','FontSize',14);
% xlabel('Expected reacation time(ms)','FontSize', 14);
% hold on
% xline(Mean_dude_random,'--r');
% hold on
% xline(Median_dude_random,'-- b');
% legend('posterior','Mean','Median');



%% TaskA3..................................................................
ZLog_Y=  log(Y);
sample_mean = zeros(1,max(ind));
for i=1:max(ind)
    sample_mean(i)=mean(ZLog_Y(ind==i));
end
figure()
scatter(1:max(ind),sample_mean,'b');
hold on;
boxplot(theta_s','Notch','off');




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



