function [tau,newcellchem,dth] = cell_gillespie(timeVec,rxnrate,state)

timeVec;
rxnrate;
state;

tau = rand(1);
newcellchem = [ones(1,8)];
dth = rand(1)*360;
dth = 360;