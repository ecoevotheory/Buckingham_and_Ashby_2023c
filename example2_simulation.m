% This code draws the evolutionary trajectories of host resistance and 
% parasite virulence for the system in Example 2

%% Define Parameters:
t_max=100;
a=10;
c=400;
alphabar=6;
sigma=10;
K=0.325;
A=1;
B=1;
cb1=0.8;
cb2=5;
resmin=0;
resmax=1;
resstart=0.5;
b=12;
alphamin=0;
alphamax=25;
alphastart=15;
res0=51;
nevol=1200;

%% First we will consider disparate host and parasite mutation rates

% Set up initial conditions:
strain_totalH = 1;
strain_totalP = 1;

% Initially we have equal numbers of susceptible and infected hosts:
init_pop = [0.1,0.1];

% Set up vectors to be used later:
Res = linspace(resmin,resmax,res0);
Alpha = linspace(alphamin,alphamax,res0);
initialH = find(Res>=resstart,1);
initialP = find(Alpha>=alphastart,1);
res_start = Res(initialH);
alpha_start = Alpha(initialP);
indexH_start = initialH;
indexP_start = initialP;

% Suppose that the parasite mutates 20 times faster than the host:
hostmutationprob=1/21;
% Now allow both traits to coevolve:
[~,~,~,~,~,~,~,RES,ALPHA,~,~] = example2_simulation_function(t_max,a,c,alphabar,sigma,b,K,A,B,cb1,cb2,resmin,resmax,res_start,alphamin,alphamax,alpha_start,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);

% Plot the host resistance trajectory:
clf
aa0=1;
RES0=log10(RES);
RES0(RES0<-aa0)=-aa0;
RES0=(RES0+aa0)/aa0;
subplot(2,2,1)
imagesc(RES0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time','interpreter','latex')
xlabel('Host resistance, $r$','interpreter','latex')
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(Res(1:(res0-1)/2:res0)*10000)/10000);
title('A')
pbaspect([1,2,1])

% Plot the parasite virulence trajectory:
ab0=1;
ALPHA0=log10(ALPHA);
ALPHA0(ALPHA0<-ab0)=-ab0;
ALPHA0=(ALPHA0+ab0)/ab0;
subplot(2,2,2)
imagesc(ALPHA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time','interpreter','latex')
xlabel('Parasite virulence, $\alpha$', 'interpreter','latex')
xlim([0,41])
set(gca,'xtick',[1,11,21,31,41],'xticklabel',[0,5,10,15,20]);
title('B')
pbaspect([1,2,1])

%% Now repeat for comparable mutation rates

% Set up initial conditions
strain_totalH = 1;
strain_totalP = 1;

% Initially we have equal numbers of susceptible and infected hosts:
init_pop = [0.1,0.1];

% Set up vectors to be used later:
Res = linspace(resmin,resmax,res0);
Alpha = linspace(alphamin,alphamax,res0);
initialH = find(Res>=resstart,1);
initialP = find(Alpha>=alphastart,1);
res_start = Res(initialH);
alpha_start = Alpha(initialP);
indexH_start = initialH;
indexP_start = initialP;

% Suppose that the parasite mutates at the same rate as the host:
hostmutationprob=1/2;
% Now allow both traits to coevolve:
[~,~,~,~,~,~,~,RES,ALPHA,~,~] = example2_simulation_function(t_max,a,c,alphabar,sigma,b,K,A,B,cb1,cb2,resmin,resmax,res_start,alphamin,alphamax,alpha_start,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);

% Plot the host resistance trajectory:
aa0=1;
RES0=log10(RES);
RES0(RES0<-aa0)=-aa0;
RES0=(RES0+aa0)/aa0;
subplot(2,2,3)
imagesc(RES0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time','interpreter','latex')
xlabel('Host resistance, $r$','interpreter','latex')
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(Res(1:(res0-1)/2:res0)*10000)/10000);
title('C')
pbaspect([1,2,1])

% Plot the parasite virulence trajectory:
ab0=1;
ALPHA0=log10(ALPHA);
ALPHA0(ALPHA0<-ab0)=-ab0;
ALPHA0=(ALPHA0+ab0)/ab0;
subplot(2,2,4)
imagesc(ALPHA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time','interpreter','latex')
xlabel('Parasite virulence, $\alpha$', 'interpreter','latex')
xlim([0,41])
set(gca,'xtick',[1,11,21,31,41],'xticklabel',[0,5,10,15,20]);
title('D')
pbaspect([1,2,1])
