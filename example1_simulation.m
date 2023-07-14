% This code draws the evolutionary trajectories of host resistance and 
% parasite infectivity for the system in Example 1

%% Set up parameters to use throughout

t_max=100;
beta0min=0;
beta0max=15;
p0=0;
p1=1;
h0=0;
h1=4;
a0=0.82;
a1=5.5;
alpha=1.5;
b=0.5;
q=0.01;
eps=0.3;

host_min=h0;
host_max=h1;
parasite_min=p0;
parasite_max=p1;

% Initial values of host and parasite traits:
hostHstart=3.2;
parasitePstart=0.03;

% There is one strain of host and one strain of parasite initially:
strain_totalH = 1;
strain_totalP = 1;

% Initial population distribution (SJ, SA, IJ, IA):
init_pop = [0.7,0,0.1,0];

%% First we consider the case where the parasite evolves at the same rate as the host

% Set up parameters and vectors:
rng(2)
res0=61;
nevol=1000;
hostHmutationprob=0.5;

%% Set up Initial Conditions

HostH = linspace(host_min,host_max,res0);
ParasiteP = linspace(parasite_min,parasite_max,res0);
initialH = find(HostH>=hostHstart,1);
initialP = find(ParasiteP>=parasitePstart,1);
hostH_start = HostH(initialH);
parasiteP_start = ParasiteP(initialP);
indexH_start = initialH;
indexP_start = initialP;

%% Allow both traits to evolve

[~,~,~,~,~,~,~,HOST,PARASITE,~,~] = example1_simulation_function(t_max,beta0min,beta0max,p0,p1,h0,h1,a0,a1,alpha,b,q,eps,hostH_start,parasiteP_start,hostHmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);

%% Make the Plot

% Plot the host resistance trajectory:
aa0=2;
HOST0=log10(HOST);
HOST0(HOST0<-aa0)=-aa0;
HOST0=(HOST0+aa0)/aa0;
subplot(2,2,3)
imagesc(HOST0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time')
xlabel('Host trait')
set(gca,'xtick',[1,31,61],'xticklabel',[0,2,4]);
title('C')
pbaspect([1,2,1])
xlim([0,res0])

% Plot the parasite virulence trajectory:
ab0=1;
PARASITE0=log10(PARASITE);
PARASITE0(PARASITE0<-ab0)=-ab0;
PARASITE0=(PARASITE0+ab0)/ab0;
subplot(2,2,4)
imagesc(PARASITE0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time')
xlabel('Parasite trait')
set(gca,'xtick',[1,31,61],'xticklabel',[0,0.5,1]);
title('D')
pbaspect([1,2,1])
xlim([1,res0])

%% Now we consider the case where the parasite evolves much faster than the host

% Set up parameters and vectors:
rng(8)
res0=201;
nevol=5000;
hostHmutationprob=1/101;

%% Set up Initial Conditions

HostH = linspace(host_min,host_max,res0);
ParasiteP = linspace(parasite_min,parasite_max,res0);
initialH = find(HostH>=hostHstart,1);
initialP = find(ParasiteP>=parasitePstart,1);
hostH_start = HostH(initialH);
parasiteP_start = ParasiteP(initialP);
indexH_start = initialH;
indexP_start = initialP;

%% Allow both traits to evolve

[~,~,~,~,~,~,~,HOST,PARASITE,~,~] = example1_simulation_function(t_max,beta0min,beta0max,p0,p1,h0,h1,a0,a1,alpha,b,q,eps,hostH_start,parasiteP_start,hostHmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);

%% Make the Plot

% Plot the host resistance trajectory:
aa0=2;
HOST0=log10(HOST);
HOST0(HOST0<-aa0)=-aa0;
HOST0=(HOST0+aa0)/aa0;
subplot(2,2,1)
imagesc(HOST0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time')
xlabel('Host trait')
set(gca,'xtick',[151,176,201],'xticklabel',[3,3.5,4]);
title('A')
pbaspect([1,2,1])
xlim([res0-61,res0])

% Plot the parasite virulence trajectory:
ab0=1;
PARASITE0=log10(PARASITE);
PARASITE0(PARASITE0<-ab0)=-ab0;
PARASITE0=(PARASITE0+ab0)/ab0;
subplot(2,2,2)
imagesc(PARASITE0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time')
xlabel('Parasite trait')
set(gca,'xtick',[1,21,41,61],'xticklabel',[0,0.1,0.2,0.3]);
title('B')
pbaspect([1,2,1])
xlim([0,61])
