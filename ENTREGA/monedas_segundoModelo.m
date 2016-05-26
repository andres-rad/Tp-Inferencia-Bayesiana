%% Inferring Number of Surveys Distributed and Return Rate Simultataneously
setenv('LD_LIBRARY_PATH', '');
clear;


% JAGS 3.2.0+ seems to work, but 3.1.0 did not

%% Sampling
n = 10;
m = 3;
k = [3,4,10];

% MCMC Parameters
nchains = 3; % How Many Chains?
nburnin = 3e3; % How Many Burn-in Samples?
nsamples = 9e5;  %How Many Recorded Samples?
nthin = 3; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Variables
datastruct = struct('k',k,'m',m,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    S.Lambda(1) = 0.5;
    S.Lambda(2) = 0.5;
    S.Lambda(3) = 0.5;
    S.c = i;
    S.Tau = [0.5 , 0.5 , 0.5];
    init0(i) = S;
end



% Use JAGS to Sample
tic
fprintf( 'Running JAGS with chains serially...\n' );
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'monedas_vector.txt'), ...
    init0, ...
    'doparallel' , doparallel, ...
    'nchains', nchains,...
    'nburnin', nburnin,...
    'nsamples', nsamples, ...
    'thin', nthin, ...
    'monitorparams', {'Theta', 'Tau',  'c', 'Lambda'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0 , ...
    'workingdir' , 'tmpjags' );
toc
    
%% Analysis

%TauAn=reshape(samples.Tau,1,[],3);
%TauPosta = (cAn==1).*TauAn(1,:,1) + (cAn==2).*TauAn(1,:,2) + (cAn==3).*TauAn(1,:,3) ;
%ThetaAn=reshape(samples.Theta,1,[],3);

%Datos generales para graficar.
cAn=reshape(samples.c,1,[]);
figure(10);clf;hold on;
ylimite = [0 15];
set(gcf,'units','norm','pos',[.2 .2 .9 .5],'paperpositionmode','auto');
nbins = 100;
wbin = 1/nbins;
binCenters = wbin/2:wbin:1-wbin/2;
bins = wbin:wbin:1;
count = nchains * nsamples;

%Procesamiento de datos
for i = 1:m
    ThetaAn(1,:,i) = reshape(samples.Theta(:,:,i), 1, []); 
end 

for i = 1:m
    TauAn(1,:,i) = reshape(samples.Tau(:,:,i), 1, []); 
    TauInter(:,i) = squeeze(TauAn(1,:,i));
end

Tau1 = TauInter(:,1);
Tau1 = Tau1(cAn == 1);
Tau1 = reshape(Tau1, 1, []);

Tau2 = TauInter(:,2);
Tau2 = Tau2(cAn == 2);
Tau2 = reshape(Tau2, 1, []);


Tau3 = TauInter(:, 3);
Tau3 = Tau3(cAn == 3);
Tau3 = reshape(Tau3, 1, []);


for i = 1:m
    LambdaAn(1,:,i) = reshape(samples.Lambda(:,:,i), 1, []); 
    LambdaInter(:,i) = squeeze(LambdaAn(1,:,i));
end

Lambda1 = LambdaInter(:,1);
Lambda1 = Lambda1(cAn ~= 1);
Lambda1 = reshape(Lambda1, 1, []);

Lambda2 = LambdaInter(:,2);
Lambda2 = Lambda2(cAn ~= 2);
Lambda2 = reshape(Lambda1, 1, []);

Lambda3 = LambdaInter(:,3);
Lambda3 = Lambda3(cAn ~= 3);
Lambda3 = reshape(Lambda3, 1, []);


%Ploteo de las densidades de Theta para cada moneda.

%Moneda 1
subplot(131);hold on; h1 = gca;
set(h1, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h1_sinN = histc(ThetaAn(1,:,1), bins);
prob1 = h1_sinN / (count * wbin);
bar(binCenters, prob1, 'hist');
title('Theta moneda 1', 'fontsize', 16);
ylim(ylimite);
xlabel('Theta');
ylabel('Densidad');

%Superpongo gr�fico de Lambda 1 

%subplot(337);hold on; h7 = gca;
%set(h7, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h7_sinN = histc(Lambda1, bins);
prob7 = h7_sinN / (size(Lambda1, 2) *wbin);
plot(binCenters, prob7,'r', 'linewidth', 1.4);
%legend('Lambda_1');
%title('Lambda moneda 1', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');


%Superpongo gr�fico de Tau 1

%subplot(334);hold on; h4 = gca;
%set(h4, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h4_sinN = histc(Tau1, bins);
prob4 = h4_sinN / (size(Tau1, 2) * wbin);
plot(binCenters, prob4,'g--', 'linewidth', 1.4);
%legend('\tau_1');
%title('Tau moneda 1', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');

legend('\theta_1', '\lambda_1', '\tau_1');


%Moneda 2
subplot(132);hold on; h2 = gca;
set(h2, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h2_sinN = histc(ThetaAn(1,:,2), bins);
prob2 = h2_sinN / (count * wbin);
bar(binCenters, prob2, 'hist');
title('Theta moneda 2', 'fontsize', 16);
ylim(ylimite);
xlabel('Theta');
ylabel('Densidad');

%Superpongo gr�fico Lambda 2

%subplot(338);hold on; h8 = gca;
%set(h8, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h8_sinN = histc(Lambda2, bins);
prob8 = h8_sinN / (size(Lambda2, 2) *wbin);
plot(binCenters, prob8,'r', 'linewidth', 1.4);
%title('Lambda moneda 2', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');

%Superpongo gr�fico Tau 2

%subplot(335);hold on; h5 = gca;
%set(h5, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h5_sinN = histc(Tau2, bins);
prob5 = h5_sinN / (size(Tau2, 2) *wbin);
plot(binCenters, prob5, 'g--', 'linewidth', 1.4);
%title('Tau moneda 2', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');

legend('\theta_2', '\lambda_2', '\tau_2');

%Moneda 3

subplot(133);hold on; h3 = gca;
set(h3, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h3_sinN = histc(ThetaAn(1,:,3), bins);
prob3 = h3_sinN / (count *wbin);
bar(binCenters, prob3, 'hist');
title('Theta moneda 3', 'fontsize', 16);
ylim(ylimite);
xlabel('Theta');
ylabel('Densidad');


%Superpongo gr�fico Lambda 3

h9_sinN = histc(Lambda3, bins);
prob9 = h9_sinN / (size(Lambda3, 2) *wbin);
plot(binCenters, prob9, 'r', 'linewidth', 1.4);
%title('Lambda moneda 3', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');

%Superpongo gr�fico Tau 3

%subplot(336);hold on; h6 = gca;
%set(h6, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h6_sinN = histc(Tau3, bins);
prob6 = h6_sinN / (size(Tau3, 2) *wbin);
plot(binCenters, prob6, 'g--', 'linewidth', 1.4);
%title('Tau moneda 3', 'fontsize', 16);
%ylim(ylimite);
%xlabel('Tau');
%ylabel('Densidad');

legend('\theta_3', '\lambda_3', '\tau_1');
%Ploteo de la distribucion de c.

c_total = reshape(samples.c, 1, []);
figure(13);clf; hold on;
title('Densidad variable categorica', 'fontsize', 16);
axis square;
hist(c_total);


stats2=stats;
stats2.mean.Tau = [mean(Tau1) mean(Tau2) mean(Tau3)];
stats2.mean.Lambda = [mean(Lambda1) mean(Lambda2) mean(Lambda3)];
stats2.std.Lambda = [std(Lambda1) std(Lambda2) std(Lambda3)];
stats2.std.Tau = [std(Tau1) std(Tau2) std(Tau3)];





