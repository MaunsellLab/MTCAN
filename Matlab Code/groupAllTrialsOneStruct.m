clear 

fileName = 'Meetz_2022_0211_MTC.dat';

header = readLLFile('i', fileName);
trials = {};
for i = 1:header.numberOfTrials
    trials{i} = readLLFile('t', i);
end

% for t = 1:length(trials)
%     trials{1,t}.spikeData = -ones([5000, 3]);
% end

cd '../Matlab Data'/
fileName = strrep(fileName, '.dat', '.mat');
save(fileName, 'trials', 'header')
% use code below to save as matlab 7.3 version
% save(fileName, 'trials', 'header', '-v7.3');
