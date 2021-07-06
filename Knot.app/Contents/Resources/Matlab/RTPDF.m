function RTPDF(dParams, subplotIndex, file, trials, indices, trialStructs)
% Plot RT PDF for task plugins.  This version expect indices to be logical arrays
%
  subplot(dParams.plotLayout{:}, subplotIndex);
  cla;
  correctRTs = [trials(indices.correct).reactTimeMS];
  earlyRTs = [trials(indices.fa).reactTimeMS];
  missRTs = [trials(indices.miss).reactTimeMS];
  allRTs = [correctRTs, earlyRTs, missRTs];
  if ~isempty(allRTs)
    cdfplot(allRTs);
  end
  % Display behavioral d' and criterion.  For this analysis we simply use the reaction time window, not 
  % fitting a response window
  earlyTimeMS = sum([trialStructs(indices.fa).preStimMS] + earlyRTs);
  earlyTimeMS = earlyTimeMS + sum([trialStructs(indices.correct).preStimMS] + file.tooFastMS);
  earlyTimeMS = earlyTimeMS + sum([trialStructs(indices.miss).preStimMS] + file.reactMS);
  releaseRateS = sum(indices.fa) / earlyTimeMS * 1000.0;
  pFA = 1.0 - exp(-releaseRateS * (file.reactMS - file.tooFastMS) / 1000.0);
  numHits = sum(indices.correct);
  pH = numHits / (numHits + sum(indices.miss));
  if (pFA > 0 && pFA < 1 && pH > 0 && pH < 1)
    [dP, ~, beta] = dPrime(pH, pFA);
    betaStr = '\beta';
    text(0.05, 0.95, sprintf('%.2f spont/s\nd'': %.1f\n%s:  %.1f', releaseRateS, dP, betaStr, beta), ...
      'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
  end
  % set the scale and add time markers for the stimulus and response limits
  if isfield(file, 'reactMS')
    timeLimit = file.reactMS;
  elseif isfield(file, 'responseLimitMS')
    timeLimit = file.responseLimitMS;
  else
    timeLimit = 1000;
  end
  set(gca, 'XLim', [-1000 timeLimit], 'YLim', [0 1]);
  hold on;
  yLimits = get(gca, 'YLim');
  plot([0 0], yLimits, 'k');
  if isfield(file, 'tooFastMS')
    plot(double(file.tooFastMS) * [1 1], yLimits, '--', 'Color', 0.5 * [0 1 0]);
  end
  if isfield(file, 'rewardedLimitMS')
   plot(double(file.rewardedLimitMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
  elseif isfield(file, 'reactMS')
     plot(double(file.reactMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
  end
  if isfield(file, 'responseLimitMS')
    plot(double(file.responseLimitMS) * [1 1], yLimits, '--', 'Color', 0.5 * [1 0 0]);
  end
  if isfield(file, 'kernelRTMinMS')
    plot(double(file.kernelRTMinMS) * [1 1], yLimits, ':', 'Color', 0.5 * [0 1 0]);
  end
  if isfield(file, 'kernelRTMaxMS')
    plot(double(file.kernelRTMaxMS) * [1 1], yLimits, ':', 'Color', 0.5 * [1 0 0]);
  end
  title('Cumulative Reaction Times');
  xlabel('');
  ylabel('');
end
