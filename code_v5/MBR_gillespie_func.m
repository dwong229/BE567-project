function [timeVec,MBRstate]= MBR_gillespie_func(rxnrate,init,simIterations,attractant)

%Initiate cells on MBR

MBRstate = struct('posn',[0 0 306*rand(1)],'cellposn',[]); %in fixed/world frame
numcell = 4;

%initialize variables
timeVec(1) = 0;
rxncell(1) = 0;


% MBR geometry:sq of sizes 40:
%place cells at center of each edge
MBRstate.cellposn(1:4,1:2) = [20,0;0,-20;-20,0;0,20]; %in MBR frame
%MBRstate.cellposn(1:4,3) = [360*rand(4,1)];
MBRstate.cellposn(1:4,3) = [0 -90 -180 -270];
    


% initial chem state:
cellstate = repmat(init,[1,1,numcell]); %cellstate(timeidx,chem,cell)

for i = 2:simIterations

    % generate a rand no to determin which cell reacts
    rxncell(i)= ceil(numcell*rand(1));

    %determine what reaction occurs and return new state and change in cell
   %angle and tau time
    [tau,newcellchem,dth] = cell_gillespie(timeVec,rxnrate,cellstate(end,:,rxncell(i)));
    timeVec(i) = timeVec(i-1) + tau;

    %update all cell states
    cellstate(end+1,:,:) = cellstate(end,:,:);
    cellstate(end+1,:,rxncell(i)) = newcellchem;
         
    %update MBR state
    % viscosity constants
    kt = 1; % p/kt
    kr = 1;
    
    bx = MBRstate.cellposn(:,1);
    by = MBRstate.cellposn(:,2);
    th = MBRstate.cellposn(:,3);
    
    dxdt_f = kt*sum(cosd(th)); %mbr frame
    dydt_f = kt*sum(sind(th)); %mbr frame
    dadt_f = kr*sum(bx.*sind(th) + by.*cosd(th));
    
    dxdt = (dxdt_f)*cosd(MBRstate.posn(i-1,3)) - dydt_f*sind(MBRstate.posn(i-1,3));
    dydt = (dxdt_f)*sind(MBRstate.posn(i-1,3)) + dydt_f*cosd(MBRstate.posn(i-1,3));
    dadt = dadt_f;
    
    % update position of MBR
    MBRstate.posn(i,:) = MBRstate.posn(i-1,:) + tau*[dxdt dydt dadt];
    MBRstate.posn(i,3) = mod(MBRstate.posn(i,3),360);
    
    % update tumbled bacterium
    MBRstate.cellposn(:,3) = MBRstate.cellposn(:,3);
    MBRstate.cellposn(rxncell(i),3) = mod(MBRstate.cellposn(rxncell(i),3) + dth,360);
    
    % dermine if ligand sensed more stuff
    
    
end
