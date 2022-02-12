numChannels = 2;
for t = 1:length(trials)
    trials{1,t}.spikeData = struct();
    for i = 1:numChannels
    trials{1,t}.spikeData.(['chan' num2str(i)]) = -ones([5000, 2]);
    end
end


fileName = 'Meetz_2022_0114_MTNAN3.dat';

header = readLLFile('i', fileName);
trials = struct();
for i = 1:header.numberOfTrials
    trials.trials{i} = readLLFile('t', i);
end

cd '../Matlab Data'/
fileName = strrep(fileName, '.dat', '.mat');
save(fileName, 'trials', 'header');


for t = 1:length(trials)
    trials{1,t}.spikeData = struct();
    for i = 1:trials{1,t}.trial.data.targetOnTimeMS
    trials{1,t}.spikeData(i).channel = -1;
    trials{1,t}.spikeData(i).unit = -1;
    trials{1,t}.spikeData(i).timeStamp = -1;
    end
end

for i = 1:10
    eventCode = NEV(i,2);
    index = find(eventDict.codes == eventCode);
    disp(eventDict.names(index))
end
