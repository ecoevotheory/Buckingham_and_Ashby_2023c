function [hostH_end,parasiteP_end,end_pop,strain_totalH,strain_totalP,indexH_end,indexP_end,HOST,PARASITE,DISPREV,NVEC] = example1_simulation_function(t_max,beta0min,beta0max,p0,p1,h0,h1,a0,a1,alpha,b,q,eps,hostH_start,parasiteP_start,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol)

% This function runs an evolutionary simulation for a fixed number of
% timesteps.

parasite_min=p0;
parasite_max=p1;
host_min=h0;
host_max=h1;

eqtol = 1e-3;
exttol = 1e-5;

% Set up vectors to be used later:
HostH = linspace(host_min,host_max,res0);
ParasiteP = linspace(parasite_min,parasite_max,res0);
HOST = zeros(nevol,res0);
PARASITE = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
NVEC = zeros(nevol,1);

% Initial conditions:
hostH_current = hostH_start;
parasiteP_current = parasiteP_start;
indexH_current = indexH_start;
indexP_current = indexP_start;

% Each value of ievol is one evolutionary timestep.
for ievol=1:nevol
    disp(ievol)
    
    % Find the ecological equilibrium:
    [~,SJ1,~,IJ1,~,~] = example1_eco_dynamics_function(t_max,beta0max,beta0min,p0,p1,h0,h1,a0,a1,alpha,b,q,eps,hostH_current,parasiteP_current,eqtol,init_pop,strain_totalH,strain_totalP);
    
    % Re-format this output into a single matrix:
    SJ=zeros(strain_totalH,strain_totalP);
    IJ=zeros(strain_totalH,strain_totalP);
    for j=1:strain_totalH
        for k=1:strain_totalP
            SJ(j,k) = SJ1(end,j+(k-1)*strain_totalH);
            IJ(j,k) = IJ1(end,j+(k-1)*strain_totalH);
        end
    end
    Nhost = SJ+IJ;
    
    % Remove extinct classes
    Nhostrows=sum(Nhost,2);
    Nhosttotal=sum(Nhost,'all');
    
    % See if any host strains go extinct:
    extinct = (Nhostrows/Nhosttotal)<exttol;
    strain_totalH = strain_totalH-sum(extinct);
    SJ(extinct,:) = [];
    IJ(extinct,:) = [];
    indexH_current(extinct) = [];
    hostH_current(extinct) = [];
    
    % Update N:
    Nparasite=IJ;
    Nparasitecolumns=sum(Nparasite,1);
    Nparasitetotal=sum(Nparasite,'all');
    
    % See if any parasite strains go extinct:
    if Nparasitetotal<exttol
        disp("Parasite is extinct")
        strain_totalP=0;
        hostH_end=hostH_current;
        parasiteP_end=parasiteP_current;
        end_pop=init_pop;
        indexH_end=indexH_current;
        indexP_end=indexP_current;
        return
    end
    extinct1 = (Nparasitecolumns/Nparasitetotal)<exttol;
    strain_totalP = strain_totalP-sum(extinct1);
    if extinct1(1)==1 && ismember(0,extinct1)
        SJ(:,find(extinct1==0,1))=SJ(:,1);
    end
    SJ(:,extinct1) = [];
    IJ(:,extinct1) = [];
    Nparasite(:,extinct1) = [];
    indexP_current(extinct1) = [];
    parasiteP_current(extinct1) = [];
    
    % Update tracker
    Nhost=SJ+IJ;
    Nhostrows=sum(Nhost,2);
    Nparasitecolumns=sum(Nparasite,1);
    Nhosttotal=sum(Nhost,'all'); 
    Nparasitetotal=sum(Nparasite,'all'); 
    
    % Proportion of hosts of each strain
    HOST(ievol,indexH_current) = Nhostrows./Nhosttotal;
    % Proportion of parasites of each strain
    PARASITE(ievol,indexP_current) = Nparasitecolumns./Nparasitetotal;
    % Proportion of individuals who have the disease
    DISPREV(ievol) = (sum(IJ,'all'))/Nhosttotal;
    % Total population density:
    NVEC(ievol) = Nhosttotal;
    
    Nhostvector=zeros(strain_totalH*strain_totalP,1);
    Nparasitevector=zeros(strain_totalH*strain_totalP,1);
    for j=1:strain_totalH
        for k=1:strain_totalP
            Nhostvector(j+(k-1)*strain_totalH)=Nhost(j,k);
            Nparasitevector(j+(k-1)*strain_totalH)=Nparasite(j,k);
        end
    end
    
    % If the mutation occurs in the host trait:
    if (rand<hostmutationprob)
        weightedprob = Nhostvector/sum(Nhostvector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locH=0;
        for j=1:strain_totalH
            for k=1:strain_totalP
                if mutator_loc==j+(k-1)*strain_totalH
                   mutator_locH=j;
                end
            end
        end
        mutator = indexH_current(mutator_locH);
 
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexH_current)) % New strain
            strain_totalH = strain_totalH+1;
            % Update vector of host trait values:
            hostH_current_again=NaN(1,length(hostH_current)+1);
            for i=1:length(hostH_current)
                hostH_current_again(1,i)=hostH_current(i);
            end
            hostH_current_again(1,end)=HostH(mutant);
            hostH_current=hostH_current_again;
            % Update vector of host trait value indices:
            indexH_current_again=NaN(1,length(indexH_current)+1);
            for i=1:length(indexH_current)
                indexH_current_again(1,i)=indexH_current(i);
            end
            indexH_current_again(1,end)=mutant;
            indexH_current=indexH_current_again;
            % Add a small population with the new trait value:
            SJ_again=NaN(size(SJ,1)+1,size(SJ,2));
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(end,:)=SJ(mutator_locH,:)/10;
            SJ=SJ_again;
           
            IJ_again=NaN(size(IJ,1)+1,size(IJ,2));
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(end,:)=IJ(mutator_locH,:)/10;
            IJ=IJ_again;
           
        end
    
    
    % If the mutation occurs in the parasite trait:
    else 
        weightedprob = Nparasitevector/sum(Nparasitevector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locP=0;
        for j=1:strain_totalH
            for k=1:strain_totalP
                if mutator_loc==j+(k-1)*strain_totalH
                   mutator_locP=k;
                end
            end
        end
        mutator = indexP_current(mutator_locP);
    
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexP_current)) % New strain
            strain_totalP = strain_totalP+1;
            % Update vector of parasite trait values:
            parasiteP_current_again=NaN(1,length(parasiteP_current)+1);
            for i=1:length(parasiteP_current)
                parasiteP_current_again(1,i)=parasiteP_current(i);
            end
            parasiteP_current_again(1,end)=ParasiteP(mutant);
            parasiteP_current=parasiteP_current_again;
            % Update vectors of parasite trait value indices:
            indexP_current_again=NaN(1,length(indexP_current)+1);
            for i=1:length(indexP_current)
                indexP_current_again(1,i)=indexP_current(i);
            end
            indexP_current_again(1,end)=mutant;
            indexP_current=indexP_current_again;
            % Add small populations with the new trait value to the
            % population:
            SJ_again=NaN(size(SJ,1),size(SJ,2)+1);
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(:,end)=zeros(size(SJ,1),1);
            SJ=SJ_again;
          
            IJ_again=NaN(size(IJ,1),size(IJ,2)+1);
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(:,end)=IJ(:,mutator_locP)/10;
            IJ=IJ_again;
         
            
        end
    end
    
    % Update initial conditions
    init_pop = zeros(1,4*strain_totalH*strain_totalP);
    for i=1:strain_totalH*strain_totalP
        init_pop(4*i-1)=IJ(i);
        init_pop(4*i-3)=SJ(i);
    end
end

% Create outputs:
hostH_end=hostH_current;
parasiteP_end=parasiteP_current;
end_pop=init_pop;
indexH_end=indexH_current;
indexP_end=indexP_current;

end