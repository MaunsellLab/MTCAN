% function [] = plotOneChannel(channel, file, trials)
% 
%   numDir = file.mapSettings.data.directionDeg.n;
%   stimDurMS = file.mapStimDurationMS.data;
%   maxNumStim = 1000;
%   histPrePostMS = 50;
%   numStim = zeros(1, numDir);
%   spikeCounts = zeros(numDir, maxNumStim);
%   spikeHists = zeros(numDir, stimDurMS + 2 * histPrePostMS);
%   for t = 1:length(trials)
%     if channel == 0
%       numMapStim = trials(t).numMap0Stim;
%       spikes = trials(t).spike0;
%     else
%       numMapStim = trials(t).numMap1Stim;
%       spikes = trials(t).spike1;
%     end
%     for s = 1:numMapStim
%       if channel == 0
%         stimDesc = trials(t).map0StimDesc(s);
%       else
%         stimDesc = trials(t).map1StimDesc(s);
%       end
%       dIndex = stimDesc.directionIndex + 1;
%       numStim(dIndex) = numStim(dIndex) + 1;
%       countStartMS = floor(trials(t).photodiodeTime + stimDesc.frameRendered * 1000.0 / stimDesc.frameRateHz);
%       countEndMS = countStartMS + stimDurMS;
%       numSpikes = sum(spikes >= countStartMS & spikes < countEndMS);
%       spikeCounts(dIndex, numStim(dIndex)) = numSpikes;
%       histStartMS = countStartMS - histPrePostMS;
%       histEndMS = countEndMS + histPrePostMS;
%       histHitBins = spikes >= histStartMS & spikes < histEndMS;
%       histIncBins = spikes(histHitBins) - histStartMS + 1;
%       spikeHists(dIndex, histIncBins) = spikeHists(dIndex, histIncBins) + 1;
%     end
%   end
%     
% 	fH = figure(1);
%   set(fH, 'units', 'inches', 'position', [26.5, 7, 8.5, 11]);
%   clf;
%   
%   % header text
%   axisHandle =   subplot(6, 4, [1, 2, 5, 6, 9, 10]);
%   set(axisHandle, 'visible', 'off');
% 	text(0.00, 1.00, 'GaborRFMap', 'FontWeight', 'bold', 'FontSize', 16, ...
%     	'horizontalAlignment', 'left',  'verticalAlignment', 'top');
% 	text(0.00, 0.90, sprintf('Subject: %s\n%s\n\nChannel %d', file.fileName, file.date, channel), 'FontSize', 12, ...
%       'horizontalAlignment', 'left',  'verticalAlignment', 'top');
% 
%   % polar plot of direction tuning
%   for d = numDir:-1:1
%     spikeMean(d) = mean(spikeCounts(d, 1:numStim(d))) * 1000.0 / stimDurMS;
%     spikeSD(d) = std(spikeCounts(d, 1:numStim(d))) * 1000.0 / stimDurMS / sqrt(numStim(d));
%   end
%   subplot(6, 4, [3, 4, 7, 8, 11, 12]);
%   polarWithErrorBars(0 :  2 * pi / numDir : 2 * pi, [spikeMean, spikeMean(1)], [spikeSD, spikeSD(1)]);
%   
%   % spike histograms
%   for d = 1:numDir
%     plotOneHist(d, spikeHists(d, :) * 1000.0 / numStim(d), stimDurMS, histPrePostMS)
%   end
%   sameAxisScaling('y', 6, 4, 13:24);
% end

% function [] = plotOneHist(d, spikeHist, stimDurMS, histPrePostMS)
% 
%   subplot(6, 4, 12 + d);
%   h = plot(smooth(spikeHist, min(0.1, 3000 / sum(spikeHist))));
%   xticks([histPrePostMS, histPrePostMS + stimDurMS])
%   xticklabels({'0', sprintf('%d', stimDurMS)});
%   title([sprintf('%d', (d - 1) * 30), char(176)]); 
%   if d == 9
%     xlabel('Time (ms)');
%     ylabel('Rate (spikes/s)');
%   end
% end

function convertAllDatsToMat()

addpath('../Matlab Code');

% Function will convert unconverted .dat files in 'Lablib Data' folder 
% to .mat file and save it in 'Matlab Data' folder

myDir = uigetdir; % sets search dir
myFiles = dir(fullfile(myDir, '*.dat')); 

% scans through files in LLDir, and converts
% them if they're not converted
for k = 1:length(myFiles)
%   fileName = erase(myFiles(k).name, '.dat');
%   fileNameMat = append(fileName, '.mat');
  fileNameMat = strrep(myFiles(1).name, '.dat', '.mat');
%   cd ..
%   cd '../Matlab Data'/
  if ~isfile([strrep(myDir, 'Lablib', 'Matlab'), '/', fileNameMat])
    convertGRF(myFiles(k).name);
%     cd ..
%     cd '../Matlab Code'/
  end    
end
end
