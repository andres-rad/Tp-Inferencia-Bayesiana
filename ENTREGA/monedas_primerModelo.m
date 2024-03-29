%% Inferring Number of Surveys Distributed and Return Rate Simultataneously
setenv('LD_LIBRARY_PATH', '');
clear;


% JAGS 3.2.0+ seems to work, but 3.1.0 did not

%% Data
n = 10;
m = 3;
k = [3,4,10];
%% Sampling
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
    S.c = i;
    S.Tau = 0.5;
    init0(i) = S;
end



% Use JAGS to Sample
tic
fprintf( 'Running JAGS with chains serially...\n' );
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, 'monedas.txt'), ...
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

%ThetaAn=reshape(samples.Theta,1,[],3);


figure(1);clf;hold on;

ylimite = [0 15];
set(gcf,'units','norm','pos',[.2 .2 .9 .5],'paperpositionmode','auto');

%Datos generales para graficar.
nbins = 100;
wbin = 1/nbins;
binCenters = wbin/2:wbin:1-wbin/2;
bins = 0:wbin:1-wbin;
count = nchains * nsamples;

%Procesamiento datos
cAn=reshape(samples.c,1,[]);
TauAn=reshape(samples.Tau,1,[]);

for i = 1:m
    ThetaAn(1,:,i) = reshape(samples.Theta(:,:,i), 1, []); 
end 

for i = 1:m-1
    LambdaAn(1,:,i) = reshape(samples.Lambda(:,:,i), 1, []);
    LambdaInter(:,i) = squeeze(LambdaAn(1,:,i));
end

Lambda1 = reshape(LambdaInter(:,1), 1, []);
Lambda2 = reshape(LambdaInter(:,2), 1, []);


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
ylabel('Count');

%Superpongo Tau
h4_sinN = histc(TauAn, bins);
prob4 = h4_sinN / (size(TauAn, 2) *wbin);
plot(binCenters, prob4,'r', 'linewidth', 1.4);

%Superpongo Lambdas
h5_sinN = histc(Lambda1, bins);
prob5 = h5_sinN / (size(Lambda1,2) * wbin);
plot(binCenters, prob5, 'g--', 'linewidth', 1.4);

h6_sinN = histc(Lambda2, bins);
prob6 = h6_sinN / (size(Lambda2, 2) * wbin);
plot(binCenters, prob6, 'k-.', 'linewidth', 1.4);

legend('\theta_1', '\tau', '\lambda_1', '\lambda_2');

%Moneda 2

subplot(132);hold on; h2 = gca;
set(h2, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h2_sinN = histc(ThetaAn(1,:,2), bins);
prob2 = h2_sinN / (count * wbin);
bar(binCenters, prob2, 'hist');
title('Theta moneda 2', 'fontsize', 16);
ylim(ylimite);
xlabel('Theta');
ylabel('Count');

%Superpongo Tau
plot(binCenters, prob4,'r', 'linewidth', 1.4);
%Superpongo Lambdas
plot(binCenters, prob5, 'g--', 'linewidth', 1.4);
plot(binCenters, prob6, 'k-.', 'linewidth', 1.4);

legend('\theta_1', '\tau', '\lambda_1', '\lambda_2');

%Moneda 3
subplot(133);hold on; h3 = gca;
set(h3, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
h3_sinN = histc(ThetaAn(1,:,3), bins);
prob3 = h3_sinN / (count *wbin);
bar(binCenters, prob3, 'hist');
title('Theta moneda 3', 'fontsize', 16);
ylim(ylimite);
xlabel('Theta');
ylabel('Count');

%Superpongo Tau
plot(binCenters, prob4,'r', 'linewidth', 1.4);
%Superpongo Lambdas
plot(binCenters, prob5, 'g--', 'linewidth', 1.4);
plot(binCenters, prob6, 'k-.', 'linewidth', 1.4);

legend('\theta_1', '\tau', '\lambda_1', '\lambda_2');

%Ploteo de la distribucion de c.

c_total = reshape(samples.c, 1, []);
figure(3);clf; hold on;
title('Densidad variable categorica', 'fontsize', 16);
axis square;
hist(c_total);




% %% Analysis
% n=reshape(samples.n,1,[]); %
% theta=reshape(samples.theta,1,[]);
% figure(50+dataset);clf;hold on;
% set(gcf,'units','norm','pos',[.2 .2 .45 .55],'paperpositionmode','auto');
% jointsize=.6;
% subplot(221);hold on;h10=gca;
% bins1=[0:15:nmax];
% bins2=[0:.03:1];
% axis([bins1(1) bins1(end) bins2(1) bins2(end)]);
% set(h10,'yaxislocation','right','box','on','fontsize',13);
% set(h10,'xtick',[],'ytick',[]);
% subplot(222);hold on; h11=gca;
% ylabel('Rate of Return','fontsize',16);
% axis([0 1 bins2(1) bins2(end)]);
% set(h11,'yaxislocation','right','ytick',[bins2(1):.2:bins2(end)],...
%     'box','on','xtick',[],'ticklength',[0 0],'fontsize',13);
% subplot(223);hold on; h12=gca;
% th=xlabel('Number of Surveys','fontsize',16);
% set(th,'rot',0,'hor','left');
% axis([bins1(1) bins1(end) 0 1]);
% set(h12,'xtick',[bins1(1):100:bins1(end)],'box','on',...
%     'ytick',[],'ticklength',[0 0],'fontsize',14);
% set(h10,'units','normalized','position',...
%     [.1 1-jointsize-.1 jointsize jointsize]);
% set(h11,'units','normalized','position',...
%     [jointsize+.1+.05 1-jointsize-.1 1-.25-jointsize jointsize]);
% set(h12,'units','normalized','position',...
%     [.1 .1 jointsize 1-.25-jointsize]);
% subplot(h10);hold on;
% ph=plot(samples.n,samples.theta,'k.');
% subplot(h11);hold on;
% count=hist(theta,bins2);
% count=count/max(count);
% ph2=barh(bins2,1-count);
% set(ph2,'facecolor','k','basevalue',1);
% subplot(h12);hold on;
% count=hist(n,bins1);
% count=count/max(count);
% ph=bar(bins1,count);
% set(ph,'facecolor','k');
% set(h11,'xlim',[1-1.2 1-0]);
% set(h12,'ylim',[0 1.2]);ph=get(gcf,'children');axes(ph(3));hold on;
% set(gca,'fontsize',14);
% 
% % Expectation
% ph=plot(mean(n),mean(theta),'rx');set(ph,'markersize',12,'linewidth',2);
% 
% % Maximum Likelihood
% cc=-inf;
% ind=0;
% for i=1:nsamples
%     logL=0;
%     for j=1:m
%         logL=logL+gammaln(n(i)+1)-gammaln(k(j)+1)-gammaln(n(i)-k(j)+1);
%         logL=logL+k(j)*log(theta(i))+(n(i)-k(j))*log(1-theta(i));
%     end;
%     if logL>cc
%         ind=i;
%         cc=logL;
%     end;
% end;
% ph=plot(n(ind),theta(ind),'go');set(ph,'markerfacecolor','w');