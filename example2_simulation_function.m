function [res_end,alpha_end,end_pop,strain_totalH,strain_totalP,indexH_end,indexP_end,RES,ALPHA,DISPREV,NVEC] = example2_simulation_function(t_max,a,c,alphabar,sigma,b,K,A,B,cb1,cb2,resmin,resmax,res_start,alphamin,alphamax,alpha_start,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol)

% This function runs an evolutionary simulation for a fixed number of
% timesteps.

eqtol = 1e-3;
exttol = 1e-5;

% Set up vectors to be used later:
Res = linspace(resmin,resmax,res0);
Alpha = linspace(alphamin,alphamax,res0);
RES = zeros(nevol,res0);
ALPHA = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
NVEC = zeros(nevol,1);

% Initial conditions:
res_current = res_start;
alpha_current = alpha_start;
indexH_current = indexH_start;
indexP_current = indexP_start;

% Each value of ievol is one evolutionary timestep.
for ievol=1:nevol
    disp(ievol)
    
    % Find the ecological equilibrium:
    alphapower=exp(-(alpha_current-alphabar).^2/sigma^2);
    [~,S1,I1,~] = example2_eco_dynamics_function(t_max,a,c,alphabar,sigma,b,K,A,B,cb1,cb2,res_current,alpha_current,alphapower,eqtol,init_pop,strain_totalH,strain_totalP);
        
    % Re-format this output into a single matrix:
    S=zeros(strain_totalH,strain_totalP);
    I=zeros(strain_totalH,strain_totalP);
    for j=1:strain_totalH
        for k=1:strain_totalP
            S(j,k) = S1(end,j+(k-1)*strain_totalH);
            I(j,k) = I1(end,j+(k-1)*strain_totalH);
        end
    end
    Nhost = S+I;
    
    % Remove extinct classes
    Nhostrows=sum(Nhost,2);
    Nhosttotal=sum(Nhost,'all');
    
    % See if any host strains go extinct:
    extinct = (Nhostrows/Nhosttotal)<exttol;
    strain_totalH = strain_totalH-sum(extinct);
    S(extinct,:) = [];
    I(extinct,:) = [];
    indexH_current(extinct) = [];
    res_current(extinct) = [];
    
    % Update N:
    Nparasite=I;
    Nparasitecolumns=sum(Nparasite,1);
    Nparasitetotal=sum(Nparasite,'all');
    
    % See if any parasite strains go extinct:
    if Nparasitetotal<exttol
        disp("Parasite is extinct")
        strain_totalP=0;
        res_end=res_current;
        alpha_end=alpha_current;
        end_pop=init_pop;
        indexH_end=indexH_current;
        indexP_end=indexP_current;
        return
    end
    extinct1 = (Nparasitecolumns/Nparasitetotal)<exttol;
    strain_totalP = strain_totalP-sum(extinct1);
    if extinct1(1)==1 && ismember(0,extinct1)
        S(:,find(extinct1==0,1))=S(:,1);
    end
    S(:,extinct1) = [];
    I(:,extinct1) = [];
    Nparasite(:,extinct1) = [];
    indexP_current(extinct1) = [];
    alpha_current(extinct1) = [];
    
    % Update tracker
    Nhost=S+I;
    Nhostrows=sum(Nhost,2);
    Nparasitecolumns=sum(Nparasite,1);
    Nhosttotal=sum(Nhost,'all'); 
    Nparasitetotal=sum(Nparasite,'all'); 
    
    % Proportion of hosts of each strain
    RES(ievol,indexH_current) = Nhostrows./Nhosttotal;
    % Proportion of parasites of each strain
    ALPHA(ievol,indexP_current) = Nparasitecolumns./Nparasitetotal;
    % Proportion of individuals who have the disease
    DISPREV(ievol) = (sum(I,'all'))/Nhosttotal;
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
    
    % If the mutation occurs in the host resistance:
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
            % Update vector of resistance trait values:
            res_current_again=NaN(1,length(res_current)+1);
            for i=1:length(res_current)
                res_current_again(1,i)=res_current(i);
            end
            res_current_again(1,end)=Res(mutant);
            res_current=res_current_again;
            % Update vector of resistance trait value indices:
            indexH_current_again=NaN(1,length(indexH_current)+1);
            for i=1:length(indexH_current)
                indexH_current_again(1,i)=indexH_current(i);
            end
            indexH_current_again(1,end)=mutant;
            indexH_current=indexH_current_again;
            % Add a small population with the new trait value:
            S_again=NaN(size(S,1)+1,size(S,2));
            for i=1:size(S,1)
                for j=1:size(S,2)
                    S_again(i,j)=S(i,j);
                end
            end
            S_again(end,:)=S(mutator_locH,:)/10;
            S=S_again;
            I_again=NaN(size(I,1)+1,size(I,2));
            for i=1:size(I,1)
                for j=1:size(I,2)
                    I_again(i,j)=I(i,j);
                end
            end
            I_again(end,:)=I(mutator_locH,:)/10;
            I=I_again;
            
        end
    
    % If the mutation occurs in the parasite virulence:
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
            % Update vector of virulence trait values:
            alpha_current_again=NaN(1,length(alpha_current)+1);
            for i=1:length(alpha_current)
                alpha_current_again(1,i)=alpha_current(i);
            end
            alpha_current_again(1,end)=Alpha(mutant);
            alpha_current=alpha_current_again;
            % Update vectors of virulence trait value indices:
            indexP_current_again=NaN(1,length(indexP_current)+1);
            for i=1:length(indexP_current)
                indexP_current_again(1,i)=indexP_current(i);
            end
            indexP_current_again(1,end)=mutant;
            indexP_current=indexP_current_again;
            % Add small populations with the new trait value to the
            % population:
            S_again=NaN(size(S,1),size(S,2)+1);
            for i=1:size(S,1)
                for j=1:size(S,2)
                    S_again(i,j)=S(i,j);
                end
            end
            S_again(:,end)=zeros(size(S,1),1);
            S=S_again;
            I_again=NaN(size(I,1),size(I,2)+1);
            for i=1:size(I,1)
                for j=1:size(I,2)
                    I_again(i,j)=I(i,j);
                end
            end
            I_again(:,end)=I(:,mutator_locP)/10;
            I=I_again;
            
        end
    end
    
    % Update initial conditions
    init_pop = NaN(1,2*strain_totalH*strain_totalP);
    for i=1:strain_totalH*strain_totalP
        init_pop(2*i)=I(i);
        init_pop(2*i-1)=S(i);
    end
end

% Create outputs:
res_end=res_current;
alpha_end=alpha_current;
end_pop=init_pop;
indexH_end=indexH_current;
indexP_end=indexP_current;

end