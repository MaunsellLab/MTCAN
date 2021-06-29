function [] = plotOneChannel(channel, file, trials)

  numEle = file.mapSettings.data.elevationDeg.n;
  numAzi = file.mapSettings.data.azimuthDeg.n;
  numGridPoints = numEle * numAzi;
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
  
  fH = figure(channel+1);
  set(fH, 'units', 'inches', 'position', [26.5, 7, 8.5, 11]);
  clf;
  
  % header text
  axisHandle =   subplot(6, 4, [1, 2, 5, 6, 9, 10]);
  set(axisHandle, 'visible', 'off');
	text(0.00, 1.00, 'GaborRF Heat Map', 'FontWeight', 'bold', 'FontSize', 16, ...
    	'horizontalAlignment', 'left',  'verticalAlignment', 'top');
	text(0.00, 0.90, sprintf('Subject: %s\n%s\n\nChannel %d', file.fileName, file.date, channel), 'FontSize', 12, ...
      'horizontalAlignment', 'left',  'verticalAlignment', 'top');
  
  %% work on this 
  % cartesian plot of speed tuning
  for i = 1:numGridPoints
    spikeMean(i) = mean(spikeCounts(i, 1:numStim(i))) * 1000.0 / stimDurMS;
    spikeSD(i) = std(spikeCounts(i, 1:numStim(i))) * 1000.0 / stimDurMS;
    spikeSEM(i) = std(spikeCounts(i, 1:numStim(i))) * 1000.0 / stimDurMS / sqrt(numStim(i));
  end
  subplot(6, 4, [3, 4, 7, 8, 11, 12]);
  X_units = (2.^(0:numSpeed-1))/file.mapSettings.data.spatialFreqCPD.minValue;
  plot(X_units, spikeMean)
  hold on
%   errorbar(X_units, spikeMean, spikeSD)
%   plot(X_units, spikeMean + spikeSEM, 'linestyle', '--')
%   plot(X_units, spikeMean - spikeSEM, 'linestyle', '--')
  
  
  
  
end
  
