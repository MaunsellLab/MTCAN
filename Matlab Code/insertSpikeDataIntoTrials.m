% Extracts channel, unit, timestamp from Nev file and inserts into apt
% trial in the .mat (originally .dat) behavioral file. 

% for MTC 
% only inserts taskEvents and taskSpike into lablib file for correct trials
% baseFileName = 'testing_220214';
% fileNameNEV = [baseFileName '_NEV.mat'];
% 
% load(fileNameNEV);
%%
% index = find(contains(taskNames, 'MTN'));
% for i = 1:length(taskEvents{1,index})
%     trialNumber = taskEvents{1,index}(i).trialStart.data;
%     for j = 1:length(trials)
%         if trials{1,j}.trialStart.data == trialNumber && trials{1,j}.extendedEOT.data == 0
%             trials{1,j}.taskEvents = taskEvents{1,index}(i);
%             % isolate rows in NEV which fall within trialStart, trialEnd window
%             trialStartS = taskEvents{1,index}(i).trialStart.timeS;
%             trialEndS = taskEvents{1,index}(i).trialEnd.timeS;
%             spikeDataIndex = find(taskSpikes{1,index}(:,3) > trialStartS & taskSpikes{1,index}(:,3) < trialEndS);
%             startInd = min(spikeDataIndex);
%             endInd = max(spikeDataIndex);
%             % trials{1,j}.spikeData = taskSpikes{1,index}(startInd:endInd, :);
%             trials{1,j}.spikeData = struct();
%             trials{1,j}.spikeData.channel = taskSpikes{1,index}(startInd:endInd,1);
%             trials{1,j}.spikeData.unit = taskSpikes{1,index}(startInd:endInd,2);
%             trials{1,j}.spikeData.timeStamp = taskSpikes{1,index}(startInd:endInd,3);
%         end
%     end
% end

% % GRF
index = find(contains(taskNames, 'GRF'));
for i = 1:length(taskEvents{1,index})
    trialNumber = taskEvents{1,index}(i).trialStart.data;
    for j = 1:length(trials)
        if trials{1,j}.trialStart.data == trialNumber && trials{1,j}.trialEnd.data == 0
            trials{1,j}.taskEvents = taskEvents{1,index}(i);
            % isolate rows in NEV which fall within trialStart, trialEnd window
            trialStartS = taskEvents{1,index}(i).trialStart.timeS;
            trialEndS = taskEvents{1,index}(i).trialEnd.timeS;
            spikeDataIndex = find(taskSpikes{1,index}(:,3) > trialStartS & taskSpikes{1,index}(:,3) < trialEndS);
            startInd = min(spikeDataIndex);
            endInd = max(spikeDataIndex);
            % trials{1,j}.spikeData = taskSpikes{1,index}(startInd:endInd, :);
            trials{1,j}.spikeData = struct();
            trials{1,j}.spikeData.channel = taskSpikes{1,index}(startInd:endInd,1);
            trials{1,j}.spikeData.unit = taskSpikes{1,index}(startInd:endInd,2);
            trials{1,j}.spikeData.timeStamp = taskSpikes{1,index}(startInd:endInd,3);
        end
    end
end


% fileNameNew = strrep(fileName, '.mat', '_withNEV.mat');
% save(fileNameNew, 'trials', 'header', '-v7.3');