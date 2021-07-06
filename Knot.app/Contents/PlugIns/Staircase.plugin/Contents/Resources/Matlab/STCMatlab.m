function dParams = STCMatlab(dParams, file, trials)

% Invoked at the end of every trial.  As a free-standing function that is instantiated on each call,
% it has no way to store static variables across calls.  Instead, such values are stored in a struct, dParams.  dParams
% arrives as an empty matrix on the first call, so the first call can be identified in this way.  By returning
% dParam with essential values, they can be recovered in each call.

  if nargin < 2
    file = []; trials = [];
  end
  [inited, dParams, file] = checkInitialization(dParams, file);
  if inited
    drawText('Staircase', dParams);
    return;
  end

  % STC can use either visual or optogenetic stimulation.  If the users changes during the day, we should only
  % analyze trials of the most recent sort.  To do this, we remove any trials that are not of the same stimulus
  % type as the most recent trial

  trialStructs = [trials(:).trial];                       % extract trial structures
  validIndices = [trialStructs(:).stimulusType] == trialStructs(end).stimulusType;
  trials = trials(validIndices);
  trialStructs = trialStructs(validIndices);
  file.trials = size(trials, 2);                          % update the number of trials
  
  % all the functions need to have the current trial ends
%   indices.correct = find([trials(:).trialEnd] == 0);      % correct trials
%   indices.fa = find([trials(:).trialEnd] == 1);           % false alarm trials
%   indices.miss = find([trials(:).trialEnd] == 2);         % miss trials

  indices.correct = [trials(:).trialEnd] == 0;           % correct trials
  indices.fa = [trials(:).trialEnd] == 1;                % false alarm trials
  indices.miss = [trials(:).trialEnd] == 2;              % miss trials

  % Call the plotting functions. Each creates a separate subplot in the figure.

  drawText('Staircase', dParams, file, trials, indices);
  dParams = RTHistogram(dParams, 7, file, trials, indices);
  outcomesOverTrial(dParams, 9, file, trials, indices, trialStructs);
  RTPDF(dParams, 4, file, trials, indices, trialStructs);
  cumulativeTime(dParams, 3, trials);
  outcomesOverDay(dParams, 6, trials, indices);
  thresholdsOverDay(dParams, trials);
  plotsvStimulus(dParams, trials, trialStructs)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots versus Stimulus Intensity
%%% These include Number of trials, Hit Rates, and Mean Reaction Times
function plotsvStimulus(dParams, trials, trialStructs)

  kPlotBins = 11;
  numTrials = length(trials);
  if trials(end).trial.stimulusType == 0      % kGaborType
    maxValue = trials(end).blockStatus.maxContrastPC;
  else                                        % kOptoType
    maxValue = trials(end).blockStatus.maxPowerMW;
  end
  bins = zeros(numTrials, 1);
  for b = 1:numTrials
    bins(b) = round(trialStructs(b).stimulusValue / maxValue * kPlotBins);
  end
  hit = [trials(:).trialEnd] == 0;
  miss = [trials(:).trialEnd] == 2;
  RTMean = zeros(kPlotBins, 1);
  RTSE = zeros(kPlotBins, 1);
  n = zeros(kPlotBins, 1);
  hits = zeros(kPlotBins, 1);
  for b = 1:kPlotBins
    binTrial = (bins(:) == b)';
    hits(b) = sum(binTrial & hit);
    n(b) = hits(b) + sum(binTrial & miss);
    RTMean(b) = mean([trials(binTrial & hit).reactTimeMS]);
    RTSE(b) = std([trials(binTrial & hit).reactTimeMS]) / sqrt(n(b));
  end
  [hitRate, pci] = binofit(hits, n);
  yNeg = hitRate - pci(:, 1);
  yPos = pci(:, 2) - hitRate;

  % Plot performance versus stimulus value
  subplot(dParams.plotLayout{:}, 5);
  cla;
  hold off;
  errorbar(1:kPlotBins, hitRate, yNeg, yPos, '-s', 'markersize', 6, 'markerfacecolor', 'blue');
  hold on;
  set(gca, 'xlim', [0, kPlotBins], 'xtick', [0, kPlotBins / 2, kPlotBins], 'XTickLabel', {'0', num2str(maxValue /2 ), num2str(maxValue)});
  title('Number of Trials');
  title('Hit Rates');
  ylabel('Percent Correct');

  % Plot RT versus stimulus value
  subplot(dParams.plotLayout{:}, 8);
  cla;
  hold off;
  errorbar(1:kPlotBins, RTMean, RTSE, '-sb', 'markersize', 6, 'markerfacecolor', 'blue', 'markeredgecolor', 'blue');
  hold on;
  set(gca, 'xlim', [0, kPlotBins], 'xtick', [0, kPlotBins / 2, kPlotBins], 'XTickLabel', {'0', num2str(maxValue /2 ), num2str(maxValue)});
  title('Mean Reaction Times');
  ylabel('RT');
  if trialStructs(end).stimulusType == 0
      xlabel('Contrast (%)');
  else
      xlabel('Power (mW)');
  end
  % Plot n versus stimulus value
  subplot(dParams.plotLayout{:}, 2);
  cla;
  hold off;
  plot(1:kPlotBins, n, '-sb', 'markersize', 6, 'markerfacecolor', 'blue', 'markeredgecolor', 'blue');
  hold on;
  set(gca, 'xlim', [0, kPlotBins], 'xtick', [0, kPlotBins / 2, kPlotBins], 'XTickLabel', {'0', num2str(maxValue /2 ), num2str(maxValue)});
  title('Number of Trials');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Threshold Progress Over Trials %%%

function thresholdsOverDay(~, trials)

  % not every trial saves questResults, so we might be called when none exist
  if ~isfield(trials, 'questResults')
      return;
  end

  qResults = [trials(:).questResults];
  qResults = qResults([qResults(:).valid] == 1);
  yNeg = [qResults(:).confidenceInterval];
  yNeg(yNeg > 10) = 0;

  subplot(4, 1, 4);
  cla;
  hold off;
  errorbar(1:length(qResults), [qResults(:).threshold], yNeg, yNeg, 'b-o', 'markerfacecolor', 'blue');
  xlabel('Trials');
  ylabel('Threshold');
  set(gca, 'ylim', [0 Inf]);
  if trials(end).trial.stimulusType == 0
      unitLabel = ' %';
  else
      unitLabel = ' mW';
  end
  title(sprintf('Threshold (%.3f%s, %d trials)', qResults(end).threshold, unitLabel, qResults(end).trials));
end
