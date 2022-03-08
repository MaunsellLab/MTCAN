close all
clear all
corrVec = [];

for i = 1:100
    fRateS = 20; % FR/s
%   stimDur = 200;
    nTrials = 10;
    nRespNeurons = 2;
    
    spikes = poissrnd(fRateS/5,nTrials,nRespNeurons);
    
    % Make Correlation Matrix
    for neuronNr1 = 1:nRespNeurons
        for neuronNr2 = 1:nRespNeurons
            if neuronNr1 == neuronNr2
                R(neuronNr1, neuronNr2) = 1;
            else neuronNr1 ~= neuronNr2;
                R(neuronNr1, neuronNr2) = 0.1;
            end
        end
    end
           
    % Cholesky decomposition of the desired correlation matrix R
    L = chol(R);
    spikes = spikes * L;
    
    % corr(spikes(:,1),spikes(:,2))
    corrVec(end+1) = corr(spikes(:,1),spikes(:,2));
end

corrMean = mean(corrVec);
corrStd = std(corrVec);
corrMed = median(corrVec);
str = {sprintf('Median %f', corrMed), sprintf('Mean %f', corrMean), sprintf('Std Dev %f', corrStd)};
histogram(corrVec)
xlim([-1 1]);
text(-0.6,20, str)
title(sprintf('Correlation Distributions for %i trials, 1000 ms stim duration', nTrials));
