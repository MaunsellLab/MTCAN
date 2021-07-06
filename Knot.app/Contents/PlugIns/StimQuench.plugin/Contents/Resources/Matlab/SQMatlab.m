function dParams = SQMatlab(dParams, file, trials)

% SQMatlab is invoked at the end of every trial.  As a free-standing function that is instantiated on each call,
% it has no way to store static across calls.  Instead, such values are stored in a struct, dParams.  dParams
% arrives as an empty matrix on the first call, so the first call can be identified in this way.  By returning
% dParam with essential values, they can be recovered in each call.

  if nargin < 2
    file = []; trials = [];
  end
  [inited, dParams, file] = checkInitialization(dParams, file);
  if inited
    drawText('StimQuench', dParams);
    return;
  end
  file.trials = size(trials, 2);              % update the number of trials

%   indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
%   indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
%   indices.miss = find([trials(:).trialEnd] == 2);         % miss trials

  indices.correct = [trials(:).trialEnd] == 0;           % correct trials
  indices.fa = [trials(:).trialEnd] == 1;                % false alarm trials
  indices.miss = [trials(:).trialEnd] == 2;              % miss trials

  trialStructs = [trials(:).trial];

  % Call the plotting functions. Each creates a separate subplot in the figure.

  drawText('StimQuench', dParams, file, trials, indices);
  RTPDF(dParams, 4, file, trials, indices, trialStructs);
  dParams = RTHistogram(dParams, 7, file, trials, indices);
  outcomesOverTrial(dParams, 9, file, trials, indices, trialStructs);
  cumulativeTime(dParams, 3, trials);
  outcomesOverDay(dParams, 6, trials, indices);
  RTvStimulus(dParams, trials, trialStructs);
  psychometric(dParams, trials, trialStructs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RT versus Stimulus Intensity %%%

function RTvStimulus(dParams, trials, trialStructs)

  hit = [trials(:).trialEnd] == 0;
  probed = [trialStructs(:).optiPowerMW] ~= 0;
  stimValueSet = unique([trials(hit).contrastPC]);
  numStim = length(stimValueSet);
  RTMean = zeros(1, numStim);
  RTSE = zeros(1, numStim);
  RTMeanProbed = zeros(1, numStim);
  RTSEProbed = zeros(1, numStim);
  for s = 1:length(stimValueSet)
      stimTrial = [trials(:).contrastPC] == stimValueSet(s);
      RTMean(s) = mean([trials(stimTrial & hit & ~probed).reactTimeMS]);
      RTSE(s) = std([trials(stimTrial & hit & ~probed).reactTimeMS]) / sqrt(sum(stimTrial & hit & ~probed));
      if sum(probed) ~= 0 % doing laser trials
          RTMeanProbed(s) = mean([trials(stimTrial & hit & probed).reactTimeMS]);
          RTSEProbed(s) = std([trials(stimTrial & hit & probed).reactTimeMS]) / sqrt(sum(stimTrial & hit & probed));
      else
          continue
      end
  end
  stimValueSet(stimValueSet == 0) = 0.1;

  subplot(dParams.plotLayout{:}, 8);
  cla;
  hold off;
  errorbar(stimValueSet, RTMean, RTSE, '-sb', 'markersize', 6, 'markerfacecolor', 'blue', 'markeredgecolor', 'blue');
  set(gca,'xscale','log', 'xgrid', 'on');
  hold on;
  errorbar(stimValueSet, RTMeanProbed, RTSEProbed, '-sb', 'markersize', 6, 'markerfacecolor', 'red', 'markeredgecolor', 'red');
  set(gca, 'XTickLabel',num2str(transpose(get(gca, 'XTick'))))          % get rid of scientific notation
  yLimits = get(gca, 'YLim');
  set(gca, 'YLim', [0 yLimits(2) * 1.05]);
  set(gca, 'XLim', [0.005 100]);
  set(gca,'XTick', [0.1 1 10 100]);
  set(gca, 'XTickLabel', {'0.1', '1', '10', '100'}) % get rid of scientific notation
  title('Mean Reaction Times');
  ylabel('RT');
  xlabel('Contrast (%)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Psychometric Function %%%

function psychometric(dParams, trials, trialStructs)

  hit = [trials(:).trialEnd] == 0;
  miss = [trials(:).trialEnd] == 2;
  stimValueSet = unique([trials(hit | miss).contrastPC]);
  numStim = length(stimValueSet);
  if numStim == 0
      return;
  end
  hits = zeros(1, numStim);
  n = zeros(1, numStim);
  probed = [trialStructs(:).optiPowerMW] ~= 0;
  for s = 1:numStim                                          % for each stim value
      stimTrial = [trials(:).contrastPC] == stimValueSet(s);
      hits(s) = sum(stimTrial & hit & ~probed);
      n(s) = hits(s) + sum(stimTrial & miss & ~probed);
      hitsProbe(s) = sum(stimTrial & hit & probed);
      nProbe(s) = hitsProbe(s) + sum(stimTrial & miss & probed);
  end
  [hitRate, pci] = binofit(hits, n);
  yNeg = hitRate - pci(:, 1)';
  yPos = pci(:, 2)' - hitRate;

  [hitRateProbed, pci] = binofit(hitsProbe, nProbe);
  yNegProbed = hitRateProbed - pci(:, 1)';
  yPosProbed = pci(:, 2)' - hitRateProbed;

  % Check and See if We are presenting a true 0 contrast
  stimValueSet(stimValueSet == 0) = 0.1;

  subplot(dParams.plotLayout{:}, 5);
  cla;
  hold off;
  errorbar(stimValueSet, hitRate, yNeg, yPos, '-s', 'markersize', 6, 'markerfacecolor', 'blue');
  hold on;
  set(gca,'xscale','log', 'xgrid', 'on');
  errorbar(stimValueSet, hitRateProbed, yNegProbed, yPosProbed, '-s', 'markersize', 6, 'markerfacecolor', 'red');
  hold on;
  % xLimits = get(gca, 'XLim');
  set(gca, 'XLim', [0.005 100]);
  set(gca,'XTick', [0.1 1 10 100]);
  set(gca, 'XTickLabel', {'0.1', '1', '10', '100'}) % get rid of scientific notation
  title(sprintf('Hit Rates, 95%% CI (n = %d - %d)', min(n), max(n)));
  ylabel('Percent Correct');
end