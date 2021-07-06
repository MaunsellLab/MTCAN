function dParams = RTHistogram(dParams, subplotIndex, file, trials, indices)
% plot the reaction time histogram for task plugin data.  This version expects the indices argument to contain
% logical arrays
%
  subplot(dParams.plotLayout{:}, subplotIndex);
  cla;
  hold on;
  correctRTs = [trials(indices.correct).reactTimeMS];
  earlyRTs = [trials(indices.fa).reactTimeMS];
  allMissRTs = [trials(indices.miss).reactTimeMS];
  missRTs = allMissRTs(allMissRTs > 0);
  if isfield(file, 'reactMS')
    timeLimit = file.reactMS;
  elseif isfield(file, 'responseLimitMS')
    timeLimit = file.responseLimitMS;
  else
    timeLimit = 1000;
  end
  edges = linspace(-1000, timeLimit, dParams.RTBins);
  if ~isempty([correctRTs, earlyRTs, missRTs])
    nCorrect = histcounts(correctRTs, edges);
    nEarly = histcounts(earlyRTs, edges);
    nMiss = histcounts(missRTs, edges);
    binSize = edges(2) - edges(1);
    bH = bar(edges(1:end - 1) + binSize / 2, [nCorrect(:), nMiss(:), nEarly(:)], 'stacked');
    set(bH, 'barWidth', 1, 'lineStyle', 'none');
    set(bH(1), 'FaceColor', [0 0 0.6]);
    set(bH(2), 'FaceColor', [0.6 0 0]);
    set(bH(3), 'FaceColor', [0.6 0 0]);
    if max([nEarly, nCorrect, nMiss] > 50)      % re-bin on next plot?
     dParams.RTBins = min([dParams.RTBins * 2, 100]);
    end
    yLimits = get(gca, 'YLim');                 % vertical line at stimulus on
    plot([0 0], yLimits, 'k');
    if isfield(file, 'tooFastMS')
      llH = plot(double(file.tooFastMS) * [1 1], yLimits, 'k--');
      set(llH, 'Color', 0.5 * [0 1 0]);
    end
    if isfield(file, 'reactMS')
      llH = plot(double(file.reactMS) * [1 1], yLimits, 'r--');
      set(llH, 'Color', 0.5 * [1 0 0]);
    elseif isfield(file, 'rewardedLimitMS')
      llH = plot(double(file.rewardedLimitMS) * [1 1], yLimits, 'r--');
      set(llH, 'Color', 0.5 * [1 0 0]);
    end
    if isfield(file, 'responseLimitMS')
      llH = plot(double(file.responseLimitMS) * [1 1], yLimits, 'r--');
      set(llH, 'Color', 0.5 * [1 0 0]);
    end
  end
  set(gca, 'XLim', [-1000 timeLimit]);
  xlabel('Time Relative to Stimulus');
  ylabel('');
  title('Reaction Times');
end
  
