function dParams = OKMatlab(dParams, file, trials)

% OKMatlab is invoked at the end of every trial. As a free-standing function that is instantiated on each call,
% it has no way to store statics across calls. Instead, such values are stored in a struct, dParams.  dParams
% arrives as an empty matrix on the first call, so the first call can be identified in this way.  By returning
% dParam with essential values, they can be recovered in each call.

  if nargin < 2
    file = []; trials = [];
  end
  [inited, dParams, file] = checkInitialization(dParams, file);
  if inited
    drawText('OKernel', dParams);
    return;
  end
  if inited                 % for offline analysis
    if nargin == 1
      drawText(dParams)
      return;
    else
      doPowerVTrial = true;
    end
  else
    doPowerVTrial = false;
  end
  
  file.trials = size(trials, 2);              % update the number of trials
  for t = 1:length(trials)
      if isempty(trials(t).reactTimeMS)
        trials(t).reactTimeMS = 10000;
      end
      if length(trials(t).reactTimeMS) > 1
          trials(t).reactTimeMS = trials(t).reactTimeMS(1);
      end
      if isempty(trials(t).trialEnd)
        trials(t).trialEnd = -1;
      end
      if length(trials(t).trialEnd) > 1
        trials(t).trialEnd = trials(t).trialEnd(1);
      end
  end

  % all the functions need to have the current trial ends

  indices.correct = [trials(:).trialEnd] == 0;           % correct trials
  indices.fa = [trials(:).trialEnd] == 1;                % false alarm trials
  indices.miss = [trials(:).trialEnd] == 2;              % miss trials

  trialStructs = [trials(:).trial];                       % trial structs extracted from trials array

  drawText('OKernel', dParams, file, trials, indices);
  OKSpecificText(file, trials);
  RTPDF(dParams, 4, file, trials, indices, trialStructs);
  outcomesOverTrial(dParams, 6, file, trials, indices, trialStructs);
  dParams  = RTHistogram(dParams, 7, file, trials, indices);
  cumulativeTime(dParams, 2, trials);
  outcomesOverDay(dParams, 3, trials, indices);

  RTvStimulus(dParams, trials, trialStructs);
  psychometric(dParams, trials, trialStructs);
  if doPowerVTrial
    oKernel(dParams, file, trials, false);
    powerVTrial(dParams, file, trials);
  else
    oKernel(dParams, file, trials, true);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Text Specific to OKernel

function OKSpecificText(file, trials)
  headerText = cell(1, 1);
  if isfield(file, 'tooFastMS') && isfield(file, 'rewardedLimitMS')
    headerText{1} = sprintf('Response window: %d -- %d ms', ...
        file.tooFastMS, file.rewardedLimitMS);
  elseif isfield(file, 'tooFastMS') && isfield(file, 'reactMS')
    headerText{1} = sprintf('Response window: %d -- %d ms', ...
        file.tooFastMS, file.reactMS);
  end
  headerText{length(headerText) + 1} = sprintf('%d ms stim, %d ms ramp', trials(1).trial.visualDurMS, ...
          trials(1).trial.visualRampDurMS);
  headerText{length(headerText) + 1} = sprintf('%d ms pulses', round(trials(1).trial.pulseDurMS));
  if isfield(file, 'kernelRTMinMS') && isfield(file, 'kernelRTMaxMS')
      headerText{length(headerText) + 1} = sprintf('Kernel RT limits: %d -- %d ms', file.kernelRTMinMS, file.kernelRTMaxMS);
  end
  text(0.00, 0.3, headerText, 'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display the optogenetic kernel %%%

function oKernel(~, file, trials, plotLastOStim)
  % if we've already plotted some trials and this one isn't stimulated or usuable, do nothing;
  if figureHasSubplotTitleString(1, 'Kernel')         % if we've already plotted some kernels
      if (trials(end).trialEnd ~= 0 && trials(end).trialEnd ~= 2) || trials(end).meanPowerMW == 0
          return
      end
  end
  startTimeMS = -200;
  endTimeMS = 200;
  meanPowers = [trials(:).meanPowerMW];
  stimIndices = meanPowers > 0;
  trialEnds = [trials(:).trialEnd];
  trialProfile = getStimProfiles(trials(end), startTimeMS, endTimeMS, false);
  if plotLastOStim
    doOneKernelPlot(9, trialProfile, 'stim', startTimeMS, endTimeMS, 'Optical Stim Last Trial', 'Power (mW)', 0, 0, false);
  end
  hitIndices = stimIndices & trialEnds == 0;
  missIndices = stimIndices & trialEnds == 2;
  % If we have kernelRTMax and Min, move excluded RTs from hits to misses
  if isfield(file, 'kernelRTMinMS') && isfield(file, 'kernelRTMaxMS')
      RTs = [trials(:).reactTimeMS];
      earlyIndices = hitIndices & RTs < file.kernelRTMinMS;
      lateIndices = hitIndices & RTs >= file.kernelRTMaxMS;
      hitIndices = hitIndices & ~(earlyIndices | lateIndices);
      missIndices = missIndices | lateIndices;
  end
  % Use the hit and miss indices to construct the profiles
  if sum(hitIndices) > 0
      hitProfiles = getStimProfiles(trials(hitIndices), startTimeMS, endTimeMS, true);
      hitMean = mean(hitProfiles);
  else
      hitMean = [];
  end
  if sum(missIndices) > 0
      missProfiles = getStimProfiles(trials(missIndices), startTimeMS, endTimeMS, true);
      missMean = mean(missProfiles);
  else
      missMean = [];
  end
  % If this is a correct trial, or there is no correct plot yet, do that plot
  if trials(end).trialEnd == 0 || ~figureHasSubplotTitleString(1, 'Kernel Hit Trials')
      if ~isempty(hitMean)
          hitCI = stimCI(size(hitProfiles, 1));
          plotTitle = sprintf('Hit Kernel (n=%d)', sum(hitIndices));
          doOneKernelPlot(10, hitMean, 'stim', startTimeMS, endTimeMS, plotTitle, 'Normalized Power',  ...
                0.5 + hitCI, 0.5 - hitCI);
      end
  end
  % If this is a miss trial, or there is no miss plot yet, do that plot
  if trials(end).trialEnd == 2 || ~figureHasSubplotTitleString(1, 'Kernel Miss Trials')
      if ~isempty(missMean)
          missCI = stimCI(size(missProfiles, 1));
          plotTitle = sprintf('Miss Kernel (n=%d)', sum(missIndices));
          doOneKernelPlot(11, missMean, 'stim', startTimeMS, endTimeMS, plotTitle, 'Normalized Power', ...
                0.5 + missCI, 0.5 - missCI);
      end
  end
  if ~isempty(hitMean) && ~isempty(missMean)
      plotTitle = sprintf('Total Kernel (n=%d)', sum(hitIndices) + sum(missIndices));
      totalMean = 2.0 * (mean([hitProfiles; -missProfiles + 1]) - 0.5);
      totalCI = stimCI(size(hitProfiles, 1) + size(missProfiles, 1));
      doOneKernelPlot(12, totalMean, 'stim', startTimeMS, endTimeMS, plotTitle, 'Normalized Power', totalCI, -totalCI);
      % make sure that the hit and miss kernels (which both exist) are on the same y-axis
      subplot(4, 3, 10);
      ax = gca;
      y0 = ylim;
      subplot(4, 3, 11);
      y1 = ylim;
      theLimits = [min(y1(1), y0(1)), max(y1(2), y0(2))];
      ylim(gca, theLimits);
      ylim(ax, theLimits);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% figureHasSubplotTitleString %%%

function hasSubplot = figureHasSubplotTitleString(figureNum, title)
    
  h = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', figureNum);
  hasSubplot = false;
  for s = 1:length(h.Children)
      if ~isempty(strfind(h.Children(s).Title.String, title))
          hasSubplot = true;
          break
      end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT versus Stimulus Intensity %%%

function RTvStimulus(dParams, trials, trialStructs)

  hit = [trials(:).trialEnd] == 0;
  stimValueSet = unique([trials(hit).visualStimValue]);
  numStim = length(stimValueSet);
  noStim =  [trialStructs(:).optoIndex] == 0;
  delay0 = [trialStructs(:).optoIndex] == 1;
  delay1 = [trialStructs(:).optoIndex] == 2;
  noStimRTMean = zeros(1, numStim);
  noStimRTSE = zeros(1, numStim);
  delay0RTMean = zeros(1, numStim);
  delay0RTSE = zeros(1, numStim);
  delay1RTMean = zeros(1, numStim);
  delay1RTSE = zeros(1, numStim);
  for s = 1:length(stimValueSet)
      stimTrial = [trials(:).visualStimValue] == stimValueSet(s);
      noStimRTMean(s) = mean([trials(stimTrial & hit & noStim).reactTimeMS]);
      noStimRTSE(s) = std([trials(stimTrial & hit & noStim).reactTimeMS]) / sqrt(sum(stimTrial & hit & noStim));
      delay0RTMean(s) = mean([trials(stimTrial & hit & delay0).reactTimeMS]);
      delay0RTSE(s) = std([trials(stimTrial & hit & delay0).reactTimeMS]) / sqrt(sum(stimTrial & hit & delay0));
      delay1RTMean(s) = mean([trials(stimTrial & hit & delay1).reactTimeMS]);
      delay1RTSE(s) = std([trials(stimTrial & hit & delay1).reactTimeMS]) / sqrt(sum(stimTrial & hit & delay1));
  end
  subplot(dParams.plotLayout{:}, 8);
  cla;
  hold off;
  errorbar(stimValueSet, delay0RTMean, delay0RTSE, '-s', 'color', [0.8500, 0.3250, 0.0980], 'markersize', 6, ...
           'markerfacecolor', [0.8500, 0.3250, 0.0980], 'markeredgecolor', [0.8500, 0.3250, 0.0980]);
  set(gca,'xscale', 'log', 'xgrid', 'on');
  hold on;
  errorbar(stimValueSet, delay1RTMean, delay1RTSE, '-s', 'color', [0.6350, 0.0780, 0.1840], 'markersize', 6, ...
           'markerfacecolor', [0.6350, 0.0780, 0.1840], 'markeredgecolor', [0.6350, 0.0780, 0.1840]);
  hold on;
  errorbar(stimValueSet, noStimRTMean, noStimRTSE, '-sb', 'markersize', 6, 'markerfacecolor', 'blue', ...
            'markeredgecolor', 'blue');
  set(gca, 'XTickLabel', num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
  yLimits = get(gca, 'YLim');
  set(gca, 'YLim', [0 yLimits(2) * 1.05]);
  title('Mean Reaction Times');
  ylabel('RT');
  if trials(end).blockStatus.visualStimType == 0
      xlabel('Contrast (%)');
  else
      xlabel('Power (mW)');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power for each trial %%

function powerVTrial(dParams, ~, trials)

  meanPowerMW = [trials(:).meanPowerMW];
  subplot(dParams.plotLayout{:}, 9);
  cla;
  title('Mean Power');
  xlabel('Trials');

  yyaxis left;
	set(gca, 'yColor', 'black');
  hold on;
  plot(meanPowerMW, 'bo');
  ylabel('Power (mW)');
  maxMW = max(meanPowerMW);
  if maxMW > 0
    set(gca, 'YLim', [0 maxMW * 1.1]);
    set(gca,'YTick', [0, maxMW]);
  end

  yyaxis right
  set(gca, 'yColor', 'black');
  correctIndices = find([trials(:).trialEnd] == 0);            % correct trials
  plot(correctIndices(2:end), diff(correctIndices), 'color', [0.8500, 0.3250, 0.0980]);
  ylabel('Trials Between Hits');
  hold off;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Psychometric Function %%%

function psychometric(dParams, trials, trialStructs)

  hit = [trials(:).trialEnd] == 0;
  miss = [trials(:).trialEnd] == 2;
  stimValueSet = unique([trials(hit | miss).visualStimValue]);
  numStim = length(stimValueSet);
  if numStim == 0
    return;
  end
  delay0 = [trialStructs(:).optoIndex] == 1;
  delay0Hits = zeros(1, numStim);
  delay0N = zeros(1, numStim);
  noStimHits = zeros(1, numStim);
  noStimN = zeros(1, numStim);
  for s = 1:numStim                                          % for each stim value
      stimTrial = [trials(:).visualStimValue] == stimValueSet(s);
      delay0Hits(s) = sum(stimTrial & hit & delay0);
      delay0N(s) = delay0Hits(s) + sum(stimTrial & miss & delay0);
      noStimHits(s) = sum(stimTrial & hit & ~delay0);
      noStimN(s) = noStimHits(s) + sum(stimTrial & miss & ~delay0);
  end
  if ~isempty(delay0Hits)
      [delay0HitRate, delay0pci] = binofit(delay0Hits, delay0N);
      delay0YNeg = delay0HitRate - delay0pci(:, 1)';
      delay0YPos = delay0pci(:, 2)' - delay0HitRate;
  end
  if ~isempty(noStimHits)
      [noStimHitRate, noStimPci] = binofit(noStimHits, noStimN);
      noStimYNeg = noStimHitRate - noStimPci(:, 1)';
      noStimYPos = noStimPci(:, 2)' - noStimHitRate;
  end
%   misses = length(miss);
%   missRate = misses / (misses + length(hit));

  subplot(dParams.plotLayout{:}, 5);
  cla;
  if exist('noStimHitRate', 'var') == 1
      hold off;
      errorbar(stimValueSet, noStimHitRate, noStimYNeg, noStimYPos, '-s', 'markersize', 6, 'markerfacecolor', 'blue');
      hold on;
  end
  if exist('delay0HitRate', 'var') == 1
      errorbar(stimValueSet, delay0HitRate, delay0YNeg, delay0YPos, '-s', 'markersize', 6, ...
           'color', [0.8500, 0.3250, 0.0980], 'markerfacecolor', [0.8500, 0.3250, 0.0980]);
      hold on;
  end
  set(gca,'xscale','log', 'xgrid', 'on');
  hold on;
  set(gca, 'ylim', [0 1]);
  xLimits = get(gca, 'XLim');
  set(gca, 'XLim', [xLimits(1), min(xLimits(2), 100)]);
  set(gca, 'XTickLabel',num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
  ylabel('Percent Correct');
  title(sprintf('Hit Rates, 95%% CI (n = %d - %d)', min([delay0N, noStimN]), max([delay0N, noStimN])));
end

