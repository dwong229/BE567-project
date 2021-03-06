%% Gillespie simulation for ecoli chemotaxis and CW and CCW motor rotation
% BE 567 Project 
% Written by: Denise Wong

function [timeVec,sysstate] = ecoli_gillespie_func(rxn_rates,initial_state,simTime,plotoption)

% Inputs:
% rxn_rates : 6x1 vector of reaction rates [lambda gamma mu delta mu_p
% delta_p]
% inital_state: 4x1 vector of the inital state of the sysem [I A M P]

%% Initialization: Reaction scheme
% States
% I: transcriptionally inactive
% A: transcriptionally active state
% M: mRNA
% P: protein

%% set defaults if no rxn rates or state is provided
if isempty(rxn_rates)
phi = 1;   % I-> A   ligand binding
delta = 1; % A -> I   reduced ligand detected
gamma = 1; % A + autophosphorelation of intracellular kinase (CheA)
alpha = 1; % phosphorelation of CheY
beta = 1;  % brings to motor (more CW rotation)
mu = 1;    % dephosphorylation
end

if isempty(initial_state)
% initialize states:
I = 0;
A = 1;
AA = 1;  %CheA
AAp = 0; %CheA-p
YY = 1;  %CheY
YYp = 0; %CheY-p
Mot = 1; %CheY-p attached to motor

x = 0;
y = 0;
theta = 0;
run = 1;
initial_state.chem = [I A AA AAp YY YYp Mot];
initial_state.dyn = [x y theta run]; %start at [x,y,theta,run]
end

if isempty(simTime)
    simTime = 100; % number of reactions to run (iterations/simulation time
end


%% Unpack inputs
% rxn rates
phi = rxn_rates(1);  % I -> A
delta = rxn_rates(2);% A -> I
alpha = rxn_rates(3);% AA -> AAp
beta = rxn_rates(4);% A + AAp + YY -> A + AA + YYp
gamma = rxn_rates(5); % YYp -> YY
mu = rxn_rates(6);   % YYp ->Mot
rho = rxn_rates(7);  % Run

% initialize states:
I = initial_state.chem(1);
A = initial_state.chem(2);
AA = initial_state.chem(3);
AAp = initial_state.chem(4);
YY = initial_state.chem(5);
YYp = initial_state.chem(6);
Mot = initial_state.chem(7);
Run = initial_state.chem(8);

x = initial_state.dyn(1);
y = initial_state.dyn(2);
th = initial_state.dyn(3);

% [I A AA AAp YY YYp Mot Run]
% each row represents the result of that rxn_i: dstate = stateChange(i,:)
stateChange = zeros(length(rxn_rates),length(initial_state.chem)); % [rxn rate x state molecules]
stateChange(1,:) = [-1 1 zeros(1,6)]; % I
stateChange(2,:) = [1 -1 zeros(1,6)]; % A
stateChange(3,:) = [0 0 -1 1 zeros(1,4)]; % AA
stateChange(4,:) = [0 0 1 -1 -1 1 0 0]; % AAp + YY
stateChange(5,:) = [0 0 0 0 1 -1 0 0]; % YYp
stateChange(6,:) = [zeros(1,6) 1 1]; % YYp
stateChange(7,:) = [zeros(1,6) -1 -1]; % Run

state(1,:) = [I A AA AAp YY YYp Mot Run];
posn(1,:) = [x y th];

timeVec(1) = 0;
rxnHistory(1) = 0;

% %% Plot
% if false
%     figure
% % Initialize Plot
% % set-up plot
% % plot states
% hold on
% hstate(1) = stairs(timeVec,state(:,1),'.-','linewidth',2,'markersize',5','color','c');
% hstate(2) = stairs(timeVec,state(:,2),'.-','linewidth',2,'markersize',5','color','g');
% hstate(3) = stairs(timeVec,state(:,3),'.-','linewidth',2,'markersize',5','color','b');
% hstate(4) = stairs(timeVec,state(:,4),'.-','linewidth',2,'markersize',5','color','r');
% hold off
% 
% legend(hstate,'I','A','M','P','Location','NorthEastOutside')
% 
% xlabel('time')
% ylabel('Number')
% 
% grid on;
% box on;
% end

%% Begin simulation
i = 1;
%while timeVec(end)<simTime
while i<simTime
    i = i+1;
    % fprintf('Iteration: %3.0f\n',i)
        
    % Unpack states
    I = state(i-1,1);
    A = state(i-1,2);
    AA = state(i-1,3);
    AAp = state(i-1,4);
    YY = state(i-1,5);
    YYp = state(i-1,6);
    Mot = state(i-1,7);
    Run = state(i-1,8);
    
    x = posn(i-1,1);
    y = posn(i-1,2);
    th = posn(i-1,3);
    % Update current reaction rates:
    rxnRateNow = rxn_rates .* [I A (AA>0) (YY>0 && AAp>0) YYp>0 YYp>0 Run>0];
    
    % Compute scaling factor to compute probability of each reaction
    rxnScale = 1/sum(rxnRateNow);
    
    rxnProb = rxnScale*rxnRateNow;
    rxnBin = cumsum(rxnProb); %generate bound of bins (if  1<r1<2, then rxn = 1)
    
    % Generate random variables: r(1): rxn, r(2): time rxn occurs
    r = rand(2,1);
    
    % Compute tau and i
       
    if max(rxnProb)<1
        %disp('More than one rxn possible')
%       disp('Inequalities')
        probrange(1,:) = [0 rxnBin(1:end-1)];
        probrange(2,:) = rxnBin;
        
        % compute when r(1) is greater then lower bound
        nextRxn = find(r(1)>probrange(1,:),true,'last');
        if r(1)>probrange(2,nextRxn)
            error('probrange error')
        end
        
%         if rxnProb(nextRxn) == 0
%             disp('rxnBin(1) == 0')
%             nextRxn = nextRxn+1;
%         elseif isempty(nextRxn)
%             disp('rxnBin(1) > r(1)')
%             nextRxn = 1;
%         end
%                 
    else
        %disp('1 rxn possible only')
        nextRxn = find(rxnProb == 1,true,'first');
    end
    %fprintf('rand(1): %3.3f\n',r(1))
    %rxnProb
    %rxnBin
    %fprintf('Reaction: %3.0f \n',nextRxn)
    
    % Update state
    state(i,:) = state(i-1,:) + stateChange(nextRxn,:);
    
    % set run state to always be at most 1:
    if state(i,8) >0
       state(i,8) = 1;
    end
    
    %state(i,:)
    
    % compute time to event
    tau = rxnScale * log(1/r(2));
    % Update time vector
    timeVec(i) = timeVec(i-1) + tau;
    % timeVec(i) % for debug
    rxnHistory(i) = nextRxn;
    
    %% update dynamics
    speed = 30; %30um/s
    if Run==0
        % tumble happens, no motion and direction changes
        thnow = normrnd(58,20);
        dx = [0 0 thnow];
        %dx = [tau*30*cos(th) tau*30*sin(th) thnow];
    else
        dx = [tau*30*cos(th) tau*30*sin(th) 0];
    end
    posn(i,:) = posn(i-1,:) + dx;
    
    
    if false
    % update plot
    %set(hstate(1),'xdata',timeVec,'ydata',state(:,1))
    %set(hstate(2),'xdata',timeVec,'ydata',state(:,2))
    set(hstate(3),'xdata',timeVec,'ydata',state(:,3))
    set(hstate(4),'xdata',timeVec,'ydata',state(:,4))
    end
end
%plot now
if plotoption
figure
%title('Gillespie Simulation')
hold on
[XX,YY] = stairs(timeVec,state(:,3));
hstate(3) = plot(XX,YY,'.-','linewidth',2,'markersize',2','color','k');
[XX,YY] = stairs(timeVec,state(:,4));
hstate(4) = plot(XX,YY,'.-','linewidth',2,'markersize',2','color','g');

hold off
legend('mRNA','Protein','Location','East')

xlabel('Time')
ylabel('Level of mRNA | Protein')
end

sysstate.chem = state;
sysstate.dyn = posn;

end

