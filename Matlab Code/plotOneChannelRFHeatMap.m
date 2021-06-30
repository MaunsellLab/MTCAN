function [] = plotOneChannelRFHeatMap(channel, file, trials)

  numEle = file.mapSettings.data.elevationDeg.n;
  numAzi = file.mapSettings.data.azimuthDeg.n;
  stimDurMS = file.mapStimDurationMS.data;
  maxNumStim = 1000;
  histPrePostMS = 50;
  numStim = zeros(numEle, numAzi);
  spikeCounts = zeros(numEle, numAzi, maxNumStim);
  spikeHists = zeros(numEle, numAzi, stimDurMS + 2 * histPrePostMS);
  
  for t = 1:length(trials)
    if channel == 0
      numMapStim = trials(t).numMap0Stim;
      spikes = trials(t).spike0;
    else
      numMapStim = trials(t).numMap1Stim;
      spikes = trials(t).spike1;
    end
    for g = 1:numMapStim
      if channel == 0
        stimDesc = trials(t).map0StimDesc(g);
      else
        stimDesc = trials(t).map1StimDesc(g);
      end
      eleIndex = stimDesc.elevationIndex + 1;
      aziIndex = stimDesc.azimuthIndex + 1;
      numStim(eleIndex, aziIndex) = numStim(eleIndex, aziIndex) + 1;  
      countStartMS = floor(trials(t).photodiodeTime + stimDesc.frameRendered * 1000.0 / stimDesc.frameRateHz);    
      countEndMS = countStartMS + stimDurMS;
      numSpikes = sum(spikes >= countStartMS & spikes < countEndMS);
      spikeCounts(eleIndex,aziIndex, numStim(eleIndex, aziIndex)) = numSpikes;
      histStartMS = countStartMS - histPrePostMS;
      histEndMS = countEndMS + histPrePostMS;
      histHitBins = spikes >= histStartMS & spikes < histEndMS;
      histIncBins = spikes(histHitBins) - histStartMS + 1;
      spikeHists(eleIndex, aziIndex, histIncBins) = spikeHists(eleIndex, aziIndex, histIncBins) + 1;
    end
  end
  
  fH = figure;
  set(fH, 'units', 'inches', 'position', [26.5, 7, 8.5, 11]);
  clf;
  
  % header text
  axisHandle =   subplot(6, 4, [1, 2, 5, 6, 9, 10]);
  set(axisHandle, 'visible', 'off');
	text(0.00, 1.00, 'GaborRF Heat Map', 'FontWeight', 'bold', 'FontSize', 16, ...
    	'horizontalAlignment', 'left',  'verticalAlignment', 'top');
	text(0.00, 0.90, sprintf('Subject: %s\n%s\n\nChannel %d', file.fileName, file.date, channel), 'FontSize', 12, ...
      'horizontalAlignment', 'left',  'verticalAlignment', 'top');
  
  % Heatmap of RF spatial pref
  spikeMean = zeros(numEle, numAzi);
  for i = 1:numEle
    for j = 1:numAzi
      spikeMean(i,j) = mean(spikeCounts(i,j, 1:numStim(i,j))) * 1000.0 / stimDurMS;
    end
  end
  subplot(6, 4, [3, 4, 7, 8, 11, 12]);
  heatmap(spikeMean);
  
  % work on this part
  % Spike Histograms of individual locations
  numLocTest = numEle * numAzi;
  count = 1;
  for i = 1:numEle
    for j = 1:numAzi
      subplot(12,6, 36+count);
      spikeHist = spikeHists(i,j,:) * 1000.0 / numStim(i,j);
      plot(smooth(spikeHist, min(0.1, 3000 / sum(spikeHist))));
      xticks([histPrePostMS, histPrePostMS + stimDurMS]);
      xticklabels({'0', sprintf('%d', stimDurMS)});
      count = count + 1;
    end
  end
  sameAxisScaling('y', 12, 6, 37:36+numLocTest);
end
   
