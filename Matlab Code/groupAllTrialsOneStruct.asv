clear all

fileName = 'Meetz_2021_08_25.dat'

header = readLLFile('i', fileName);
trials = struct();
for i = 1:header.numberOfTrials
    trials.trials{i} = readLLFile('t', i);
end

cd '../Matlab Data'/
fileName = strrep(fileName, '.dat', '.mat');
save(fileName, 'trials', 'header');
