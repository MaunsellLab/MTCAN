% Extracts channel, unit, timestamp from Nev file and inserts into apt
% trial in the .mat (originally .dat) behavioral file. 

eventDict = taskDict();
trialsCount = 1;

for i = 1:length(NEV)
    if NEV(i,1) == 0 && NEV(i,2) == 54355  %trialStart event code 
        trialStartIndex = i;
        chan0Counter = 0;
        for j = trialStartIndex+1:length(NEV)
            while chan0Counter < 1
                if NEV(j,1) == 0 && ~ismember(NEV(j,2), eventDict.codes)
                    nevTrialCounter = NEV(j,2);
                    chan0Counter = chan0Counter + 1;
                end
            end
            if NEV(j,1) == 0 && NEV(j,2) == 54341 %trialEnd event code
                for x = j+1:length(NEV)
                    while chan0Counter < 2
                        if NEV(x,1) == 0
                            trialEndIndex = x;
                            for t = 1:length(trials)
                                if trials{1,t}.trialStart.data == nevTrialCounter
                                    trials{1,t}.spikeData = NEV(trialStartIndex:trialEndIndex,:);
                                    trialsCount = trialsCount + 1;
                                end
                            end
                            chan0Counter = chan0Counter + 1;
                        end 
                    end
                end
                break
            end
        end
        % break
    end
end
fileNameNew = strrep(fileName, '.mat', '_withNEV.mat');
save(fileNameNew, 'trials', 'header', '-v7.3');