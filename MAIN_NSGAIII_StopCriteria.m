% =====================================================================  
% --------------------------------------------------------------------- 
% ---------------------   JUAN CARLOS TICONA  -------------------------
% ---------- INSTITUTO DE PESQUISAS HIDRAULICAS (IPH) UFRGS  ----------
% -------------------------- OUTUBRO DE 2023 --------------------------    
% --------------------------------------------------------------------- 
% =====================================================================
% 
t = zeros(3,10);
for arquivo = 1:10
    save ('salva.mat','arquivo','t');    
    clc;
    clear;
    close all;
    load ('salva.mat');
    disp(arquivo);
    tic;
    primeira = 0;
    segunda = 0;
    Count = 0;
    Count_Max = 10;
    F1min = 0.0;
    F2min = 0.0;
    d = [];
    numPop = [];
    GD0=10000;
    
    %% Definição do Problema
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Valores utilizados no HyMOD_5p
%     nVar = 5;    %Variáveis de decisão.
%     % Parâmetros do Tank-Model a serem calibrados.
%     % X = [Smax, b, a, kf, ks]
%     %Limite inferior de cada uma das variáveis de decisao:
%     VarMin=[10.0,0.0,0.0,0.15,0.0];   
%     %Limite superior de cada uma das variáveis de decisao:
%     VarMax=[2000,7.0,1.0,1.0,0.15];    

%     VarSize = [1 nVar];         %Size of Decision Variables Matrix
%     FO = @(x) FO_HYMOD(x);       % Function Objective    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Valores utilizados no GR4J_4p
%     nVar = 6;    %Variáveis de decisão.
%     % Parâmetros do Tank-Model a serem calibrados.
%     % X = [Smax, kf, Rmax, T, So, Fo]
%     %Limite inferior de cada uma das variáveis de decisao:
%     VarMin=[0.01,-10.0,10.0,0.5,1,1]; 
%     %Limite superior de cada uma das variáveis de decisao:
%     VarMax=[1500,5.0,500.0,4.0,100,200];  

%     VarSize = [1 nVar];         %Size of Decision Variables Matrix
%     FO = @(x) FO_GR4J(x);       % Function Objective   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Valores utilizados no GR5J_7p nVar = 7;     
%     nVar = 7;    %Variáveis de decisão.
%     % Parâmetros do Tank-Model a serem calibrados.
%     % X = [Smax, kf, Rmax, T, K, So, Fo]
%     %Limite inferior de cada uma das variáveis de decisao:
%     VarMin=[0.01,-10.0,1.0,0.5,0.001,1,1]; 
%     %Limite superior de cada uma das variáveis de decisao:
%     VarMax=[2000,5.0,500.0,4.0,1.0,100,200];   

%     VarSize = [1 nVar];         %Size of Decision Variables Matrix
%     FO = @(x) FO_GR5J(x);       % Function Objective   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Valores utilizados no IPH2_7p
%     nVar = 7;    %Variáveis de decisão.
%     % Parâmetros do Tank-Model a serem calibrados.
%     % X = [Io, Ib, h, Ksup, Ksub, Rmax, a]
%     %Limite inferior de cada uma das variáveis de decisao:
%     VarMin = [10.0,0.1,0.01,0.01,10.0,0.0,0.01];  
%     %Limite superior de cada uma das variáveis de decisao:
%     VarMax = [300.0,10.0,0.99,10.0,500.0,9.0,20.0];

%     VarSize = [1 nVar];         %Size of Decision Variables Matrix
%     FO = @(x) FO_IPH2(x);       % Function Objective  
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Valores utilizados no Tank-Model 3 hidro
%     nVar = 9;    %Variáveis de decisão.
%     % Parâmetros do Tank-Model a serem calibrados.
%     % X = [H1, H2, H3, a1, a2, a3, a4, b1, b2]
%     %Limite inferior de cada uma das variáveis de decisao:
%     VarMin=[10,10,10,0.09,0.09,0.09,0.01,0.01,0.01];   
%     %Limite superior de cada uma das variáveis de decisao:
%     VarMax=[70,45,70,0.5,0.5,0.5,0.1,0.1,0.1];    

%     VarSize = [1 nVar];         %Size of Decision Variables Matrix
%     FO = @(x) FO_TANK3hidro(x);      % Function Objective   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Valores utilizados no Tank-Model 4 hidro
    nVar = 16;    %Variáveis de decisão. 
    % Parâmetros do Tank-Model a serem calibrados.
    % X = [HI1, HI2, HI3, HI4, HA1, HA2, HB1, HC1, a1, a2, b1, c1, d1, a0, b0, c0]
    %Limite inferior de cada uma das variáveis de decisão:
    VarMin=[5,10,50,100,10,10,10,10,0.09,0.09,0.09,0.01,0.001,0.01,0.01,0.001];   
    %Limite superior de cada uma das variáveis de decisão:
    VarMax=[75,70,200,500,70,45,70,70,0.5,0.5,0.5,0.1,0.01,0.1,0.1,0.1];  

    VarSize = [1 nVar];          % Size of Decision Variables Matrix
    FO = @(x) FO_TANK4hidro(x);  % Function Objective

    % Número de funções objetivas
    nObj = numel(FO(VarMin + (VarMax - VarMin)*rand));
        
    %% NSGA-III Parâmetros

    % Generating Reference Points
    nDivision = 10;
    Zr = GenerateReferencePoints(nObj, nDivision);

    
    MaxIt=300;      % Número Máximo de Iterações    
    nPop=50;        % Tamanho da população
    
    pCrossover=0.9;                         % Porcentagem de soluções obtidas do Cruzamento
    nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)
    
    pMutation=0.5;                          % Porcentagem de soluções obtidas da Mutação
    nMutation=round(pMutation*nPop);        % Number of Mutants
    
    mu=0.02;                    % Taxa de Mutação
    
    sigma=0.05*(VarMax-VarMin);  % Tamanho do passo de mutação
    
    StopCrit = zeros(MaxIt,2);
    
    %% Iniziando Parametros

    params.nPop = nPop;
    params.Zr = Zr;
    params.nZr = size(Zr,2);
    params.zmin = [];
    params.zmax = [];
    params.smin = [];

    %% Inicialização
    
    disp('Staring NSGA-III ...');
    
    empty_individual.Position = [];
    empty_individual.Cost = [];
    empty_individual.Rank = [];
    empty_individual.DominationSet = [];
    empty_individual.DominatedCount = [];
    empty_individual.NormalizedCost = [];
    empty_individual.AssociatedRef = [];
    empty_individual.DistanceToAssociatedRef = [];

    pop=repmat(empty_individual,nPop,1);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Criação das populações para os testes
%     for i=1:nPop
%         pop(i).Position=VarMin + (VarMax - VarMin)*rand;
%         pop(i).Cost=FO(pop(i).Position);
%     end
%     
%     save(sprintf('população_goias_Tank4_%d.mat',arquivo),'pop')
    
    %%
    population = load (sprintf('população_goias_IPHII_%d.mat',arquivo));
    for i=1:nPop
        pop(i).Position=population.pop(i).Position;
        pop(i).Cost=population.pop(i).Cost;
    end
    
    % Classifique a população e execute a seleção
    [pop, F, params] = SortAndSelectPopulation(pop, params);

    %% NSGA-III Loop principal
    
    for it=1:MaxIt
        
        % Crossover
        popc=repmat(empty_individual,nCrossover/2,2);
        for k=1:nCrossover/2
            
            i1=randi([1 nPop]);
            p1=pop(i1);
            
            i2=randi([1 nPop]);
            p2=pop(i2);
            
            [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
            
            popc(k,1).Cost=FO(popc(k,1).Position);
            popc(k,2).Cost=FO(popc(k,2).Position);
            
        end
        popc=popc(:);
        
        % Mutação
        popm=repmat(empty_individual,nMutation,1);
        for k=1:nMutation
            
            i=randi([1 nPop]);
            p=pop(i);
            
            popm(k).Position=Mutate(p.Position,mu,sigma,VarMin,VarMax);
            
            popm(k).Cost=FO(popm(k).Position);
            
        end
        
        % Mesclar
        pop=[pop
            popc
            popm]; % #ok
        
        % Ordenar População e execute a seleção
        [pop, F, params] = SortAndSelectPopulation(pop, params);
        
        % Store F1
        F1=pop(F{1});
        Nsol = numel(F1);
        
        % Calculo de Criterios de Parada
       [nF GD] = StopingCriteria(F1);
        StopCrit(it,:) = [nF GD];
       
        % Criterio de parada das iterações
        Melhores=zeros(nF,nVar+3);
        for i = 1:nF
%             Melhores(i,:) = [F1(i).Position (F1(i,1).Cost(1)) (F1(i,1).Cost(2)) sqrt((F1(i,1).Cost(1)-F1min)^2 + (F1(i,1).Cost(2)-F2min)^2)];
            Melhores(i,:) = [F1(i).Position (1-F1(i,1).Cost(1)) (1-F1(i,1).Cost(2)) sqrt((F1(i,1).Cost(1)-F1min)^2 + (F1(i,1).Cost(2)-F2min)^2)];
        end
        if nF > nPop-1 && sum(isnan(Melhores(:,nVar + 1))) == 0 && sum(isnan(Melhores(:,nVar + 2))) == 0 && it > 125
            Count = Count+1;
        else
            Count=0;
        end
        %% Criterio padrao CP1
        if Count > Count_Max && primeira == 0
            disp('parou CP1');            
            %% Results
            Melhores1 = sortrows(Melhores, nVar + 3);
            xlswrite(sprintf('pareto goias NSGAIII_3CP+IPH2_T2.xlsx'),Melhores1,num2str(['CP1_' num2str([arquivo it])]))
            primeira = 1;
            %break;
            t(1,arquivo) = toc;
        end
        %% Criterio padrao CP2
        if GD < 0.15
            if abs(GD0 - GD) < 0.001 && segunda == 0;
                disp('parou CP2');            
                %% Results
                Melhores2 = sortrows(Melhores, nVar + 3);
                xlswrite(sprintf('pareto goias NSGAIII_3CP+IPH2_T2.xlsx'),Melhores2,num2str(['CP2_' num2str([arquivo it])]))
                segunda = 1;
                %break;
                t(2,arquivo) = toc;
            else
            GD0 = GD;
            end
        end
        %% Criterio padrao CP3
        if it == MaxIt
            %% Results
            disp('parou CP3');
                disp('Optimization Terminated.');
            Melhores3 = sortrows(Melhores, nVar + 3);
            xlswrite(sprintf('pareto goias NSGAIII_3CP+IPH2_T2.xlsx'),Melhores3,num2str(['CP3_' num2str([arquivo it])]))
%             xlswrite(sprintf('Criterio de parada goias NSGAIII+IPH2_T2.xlsx'),Title,num2str([arquivo it]),'A1')
            xlswrite(sprintf('Medidas goias NSGAIII+IPH2_T2.xlsx'),StopCrit,num2str([arquivo it]),'A2')
            t(3,arquivo) = toc;
            break;
        end
    end
end
save ('tempo 3C_goias_NSGAIII_IPH2_T2.mat','t');