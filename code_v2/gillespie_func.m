%% Gillespie simulation
% BE 567 Homework 3, Fall 2012
% Written by: Denise Wong

function [timeVec,state] = gillespie_func(rxn_rates,initial_state,simTime,plotoption)

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

if isempty(rxn_rates)
lambda = 1;  % I -> A
gamma = 0;   % A -> I
mu = 1;      % A -> A+M
delta = 1;   % M -> 0
mu_p = 0;    % M -> M + P
delta_p = 0; % P -> 0
end

if isempty(initial_state)
% initialize states:
I = 0;
A = 1;
M = 0;
P = 0;
end

lambda = rxn_rates(1);  % I -> A
gamma = rxn_rates(2);   % A -> I
mu = rxn_rates(3);      % A -> A+M
delta = rxn_rates(4);   % M -> 0
mu_p = rxn_rates(5);    % M -> M + P
delta_p = rxn_rates(6); % P -> 0

% initialize states:
I = initial_state(1);
A = initial_state(2);
M = initial_state(3);
P = initial_state(4);

if isempty(simTime)
simTime = 100; % number of reactions to run (iterations/simulation time
end

% [I A M P]
% each row represents the result of that rxn_i: dstate = stateChange(i,:)
stateChange = [-1 1 0 0;...
    1 -1 0 0;...
    0 0 1 0;
    0 0 -1 0;
    0 0 0 1;
    0 0 0 -1];

%rxn_rate = [lambda gamma mu delta mu_p delta_p];

state(1,:) = [I A M P];

timeVec(1) = 0;
rxnHistory(1) = 0;

if false
    figure
% Initialize Plot
% set-up plot
% plot states
hold on
hstate(1) = stairs(timeVec,state(:,1),'.-','linewidth',2,'markersize',5','color','c');
hstate(2) = stairs(timeVec,state(:,2),'.-','linewidth',2,'markersize',5','color','g');
hstate(3) = stairs(timeVec,state(:,3),'.-','linewidth',2,'markersize',5','color','b');
hstate(4) = stairs(timeVec,state(:,4),'.-','linewidth',2,'markersize',5','color','r');
hold off

legend(hstate,'I','A','M','P','Location','NorthEastOutside')

xlabel('time')
ylabel('Number')

grid on;
box on;

title('Gillespie Simulation: Denise Wong')
end
%% Begin simulation
i = 1;
%while timeVec(end)<simTime
while i<simTime
    i = i+1;
    % fprintf('Iteration: %3.0f\n',i)
    
    
    % Unpack states
    I = state(i-1,1);
    A = state(i-1,2);
    M = state(i-1,3);
    P = state(i-1,4);
    
    % Update current reaction rates:
    rxnRateNow = rxn_rates .* [I A A M M P];
    
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
    %state(i,:)
    
    % compute time to event
    tau = rxnScale * log(1/r(2));
    % Update time vector
    timeVec(i) = timeVec(i-1) + tau;
    % timeVec(i) % for debug
    rxnHistory(i) = nextRxn;
    
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

end

