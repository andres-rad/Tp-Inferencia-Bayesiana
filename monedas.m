%% Inferring Number of Surveys Distributed and Return Rate Simultataneously

clear;


% JAGS 3.2.0+ seems to work, but 3.1.0 did not

%% Data
n = 10;
m = 3;
k = [10,10,5];
%% Sampling
% MCMC Parameters
nchains = 3; % How Many Chains?
nburnin = 3e3; % How Many Burn-in Samples?
nsamples = 5e4;  %How Many Recorded Samples?
nthin = 5; % How Often is a Sample Recorded?
doparallel = 0; % Parallel Option

% Assign Matlab Variables to the Observed Variables
datastruct = struct('k',k,'m',m,'n',n);

% Initialize Unobserved Variables
for i=1:nchains
    % S.Thetao = 0.5;
    S.c = i;
    S.Thetau = 0.5;
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
    'monitorparams', {'Theta', 'Thetau',  'c', 'Thetao'}, ...
    'savejagsoutput' , 1 , ...
    'verbosity' , 1 , ...
    'cleanup' , 0 , ...
    'workingdir' , 'tmpjags' );
toc
    
%% Analysis
cAn=reshape(samples.c,1,[]);
ThetauAn=reshape(samples.Thetau,1,[]);
%ThetaAn=reshape(samples.Theta,1,[],3);
for i = 1:3
    ThetaAn(1,:,i) = reshape(samples.Theta(:,:,i), 1, []); 
end 
figure(1);clf;hold on;

ylimite = [0 12000];
set(gcf,'units','norm','pos',[.2 .2 .9 .5],'paperpositionmode','auto');


subplot(131);hold on; h1 = gca;
set(h1, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
hist(ThetaAn(1,:,1), 60);
title('Theta moneda 1');
ylim(ylimite);
xlabel('Theta');
ylabel('Count');


subplot(132);hold on; h2 = gca;

set(h2, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
hist(ThetaAn(1,:,2), 60);
title('Theta moneda 2');
ylim(ylimite);
xlabel('Theta');
ylabel('Count');


subplot(133);hold on; h3 = gca;
set(h3, 'yaxislocation', 'left', 'box', 'on', 'fontsize', 13);
hist(ThetaAn(1,:,3), 60);
title('Theta moneda 3');
ylim(ylimite);
xlabel('Theta');
ylabel('Count');

figure(2);clf;hold on;
%subplot(121);hold on;
axis square;
eps=.01;bins=[0:eps:1];binsc=[eps/2:eps:1-eps/2];
count=histc(reshape(samples.Thetau,1,[]),bins);
count=count(1:end-1);
count=count/sum(count)/eps;
ph=plot(binsc,count,'k-');
set(gca,'xlim',[0 1],'box','on','fontsize',14,'xtick',[0:.2:1],'ytick',[1:ceil(max(get(gca,'ylim')))]);
set(gca,'box','on','fontsize',14);
xlabel('Rate','fontsize',16);
ylabel('Density','fontsize',16);


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
