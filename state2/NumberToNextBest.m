function stop = NumberToNextBest(optimValues, state)

persistent bestfv bestcounter

stop = false;
switch state
    case 'init'
        % Initialize variable to record best function value.
        bestfv = []; 
        
        % Initialize counter to record number of
        % local solver runs to find next best minimum.
        bestcounter = 1; 
        
        % Create the histogram.
        bar(log(bestcounter),'tag','NumberToNextBest');
        xlabel('Number of New Best Fval Found');
        ylabel('Log Number of Local Solver Runs');
        title('Number of Local Solver Runs to Find Lower Minimum')
    case 'iter'
        % Find the axes containing the histogram.
        NumToNext = ...
          findobj(get(gca,'Children'),'Tag','NumberToNextBest');
        
        % Update the counter that records number of local
        % solver runs to find next best minimum.
        if ~isequal(optimValues.bestfval, bestfv)
            bestfv = optimValues.bestfval;
            bestcounter = [bestcounter 1];
        else
            bestcounter(end) = bestcounter(end) + 1;
        end
        
        % Update the histogram.
        set(NumToNext,'Ydata',log(bestcounter))
end
