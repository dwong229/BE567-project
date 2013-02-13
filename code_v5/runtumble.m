function [runTime,tumbleTime,rtRatio,dx] = runtumble(timeVecAll,state,timeLim)
% use with ecoli simulation
% compute runTime, tumble Time and run/tumble time ratio


% time is given in seconds, need to sample at the same time each round
timeVec = timeVecAll;
timeVec(timeVecAll>timeLim) = [];
chem = state.chem(1:length(timeVec),:);
posn = state.dyn(1:length(timeVec),:);
diffTime = timeVec(2:end) - timeVec(1:end-1);
runTime = sum(diffTime(chem(1:end-1,8)>0));
tumbleTime= timeVec(end) - runTime;
rtRatio = runTime/tumbleTime;
dx = ([posn(end,1:2) - posn(1,1:2)]);

if isinf(rtRatio)
    disp('rtRatio inf')
    keyboard
end

end
