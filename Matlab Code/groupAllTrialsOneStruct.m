header = readLLFile('i', 'testing_2021_0912.dat');

trials = struct()
for i = 1:header.numberOfTrials
    trials.trials{i} = readLLFile('t', i);
end

% trials.trials = reshape(trials.trials, [length(trials.trials),1])