% E.coli Simulation Main
clear
close all
%%%%%%%%
simMode = 4; %1: one simulation, 2: repetition
%%%%%%%%
repeatSim = 1000;
simIterations = 1000;

%% simulation parameters
phi = 0;   % I-> A   ligand binding
delta = 1; % A -> I   reduced ligand detected
alpha = 0.11; % AA -> AAp
beta = 1e8;    % AAp + YY -> AA + YYp autophosphorelation of intracellular kinase (CheA) + phosphorelation of CheY
gamma = 0.031; % YYp -> YY
%gamma = 100; % YYp -> YY
mu = 0;
%mu = 0.015;      % YYp -> motor
rho = 0.015/10;     % run rate (run to tumble)
rho = (alpha+beta+gamma+mu)/10;
mu = 10; % tumble to run
rho = 1/10; % run to tumble

rxnrate = [phi delta alpha beta gamma mu rho];

% I = 1; % no chemoattractant
% A = 0;
% AA = 1;  %CheA
% AAp = 0; %CheA-p
% YY = 1;  %CheY
% YYp = 0; %CheY-p
% Mot = 1; %CheY-p attached to motor
% Run = 1;

I = 1; % no chemoattractant
A = 0;
%AA = 6700;  %CheA
AA = 2;
AAp = 0; %CheA-p
%YY = 8200;  %CheY
YY = 2;
YYp = 0; %CheY-p
Mot = 1; %CheY-p attached to motor
Run = 1;


x = 0;
y = 0;
theta = 360*rand(1);

init.chem = [I A AA AAp YY YYp Mot Run];
init.dyn = [x y theta]; %start at [x,y,theta]

%simTime = 300;

switch simMode
    
    case 1
        % Computation for single cell:
        disp('Running Simulation Once')
        % Simulate
        
        attractant = @(x) 0;
        [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
        simTime = floor(timeVec(end)*0.90);
        
        %plot simulation
        plotStates(timeVec,state)
        
    case 2
        disp('Running Repetition of Simulation')
        
        t = CTimeleft(repeatSim);
        for rep = 1:repeatSim
            t.timeleft();
            plotnow = false;
            
            % run simulation
            attractant = @(x) 0;
            [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
            
            % compute stats
            if rep == 1
                simTime = floor(timeVec(end)*0.75);
            end
            
            [~,~,rtratio,dx] = runtumble(timeVec,state,simTime);
            
            stats.dx(rep,:) = dx;
            stats.rtratio(rep) = rtratio;
        end
        save('randomwalkstats.mat','stats')
        
        % plot statistics
        figure
        hist(stats.rtratio)
        title('Histogram of run/tumble ratio')
        xlabel('run/tumble ratio')
        ylabel('frequency')
        
        figure
        for i = 1:size(stats.dx,1)
            hold on
            line([0 stats.dx(i,1)],[0 stats.dx(i,2)])
        end
        title('Random Walk Planar Displacement')
        xlabel('x (um)')
        ylabel('y (um)')
        
        avgrtratio = mean(stats.rtratio);
        stdrtratio = std(stats.rtratio);
        fprintf('Mean run-tumble ratio: %4.4f, stdev: %4.4f \n',avgrtratio,stdrtratio)
        
        avgdx = mean(stats.dx,1);
        stddx = std(stats.dx,1);
        fprintf('Mean displacement ratio: (%4.0f,%4.0f), stdev: (%4.0f,%4.0f) \n',avgdx(1),avgdx(2),stddx(1),stddx(2))
        
    case 3
        % Computation for single cell:
        disp('Running Simulation Once with Attractant')
        % Simulate
        
        attractant = @(x) attract_exp(x);
        [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
        simTime = floor(timeVec(end)*0.75);
        
        %plot simulation
        plotStates(timeVec,state)
        
        pathhandle = gcf;
        xmin = min(state.dyn(:,1));
        xmax = max(state.dyn(:,1));
        ymin = min(state.dyn(:,2));
        ymax = max(state.dyn(:,2));
        plot_attract_exp(pathhandle,[xmin xmax ymin ymax])
    
    case 4
        disp('Running Repetition of Simulation with attractant')
        
        t = CTimeleft(repeatSim);
        for rep = 1:repeatSim
            t.timeleft();
            plotnow = false;
            
            % run simulation
            
            attractant = @(x) attract_exp(x);
            [timeVec,state] =  ecoli_gillespie_func(rxnrate,init,simIterations,attractant);
            
            % compute stats
            if rep == 1
                simTime = floor(timeVec(end)*0.75);
            end
            
            [~,~,rtratio,dx] = runtumble(timeVec,state,simTime);
            
            stats.dx(rep,:) = dx;
            stats.rtratio(rep) = rtratio;
        end
        save('attractexpstats.mat','stats')
        
        % plot statistics
        figure
        hist(stats.rtratio)
        title('Histogram of run/tumble ratio')
        xlabel('run/tumble ratio')
        ylabel('frequency')
        
        figure
        for i = 1:size(stats.dx,1)
            hold on
            line([0 stats.dx(i,1)],[0 stats.dx(i,2)])
        end
        title('Random Walk Planar Displacement')
        xlabel('x (um)')
        ylabel('y (um)')
        
        avgrtratio = mean(stats.rtratio);
        stdrtratio = std(stats.rtratio);
        fprintf('Mean run-tumble ratio: %4.4f, stdev: %4.4f \n',avgrtratio,stdrtratio)
        
        avgdx = mean(stats.dx,1);
        stddx = std(stats.dx,1);
        fprintf('Mean displacement ratio: (%4.0f,%4.0f), stdev: (%4.0f,%4.0f) \n',avgdx(1),avgdx(2),stddx(1),stddx(2))
        
    
    otherwise
        disp('Invalid simMode')
end

%%
% Compute amount of time in Run vs. Tumble
% when Run == 1: run
disp('....For the last round')
[runTime,tumbleTime,rtRatio,dx] = runtumble(timeVec,state,simTime);
fprintf('Run Time: %4.1f s, Tumble Time: %4.1f s. Run/Tumble: %4.4f \n',runTime, tumbleTime,rtRatio)
fprintf('Displacement: (%4.0f,%4.0f) um \n',dx(1),dx(2))

tumble = state.chem(1:end-1,8) - state.chem(2:end,8);
tFreq = sum(tumble==1)/timeVec(end);
%rFreq = sum(tumble ==-1)/timeVec(end);
fprintf('Tumble Freq: %4.2f /s \n',tFreq)