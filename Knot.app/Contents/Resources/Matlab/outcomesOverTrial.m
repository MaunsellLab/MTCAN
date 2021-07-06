%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trial outcomes over the time course of individual trials %%%

function outcomesOverTrial(dParams, plotIndex, file, trials, indices, trialStructs)

  preStimMS = [trialStructs(:).preStimMS];
  hitTimes = preStimMS(indices.correct);
  hitRTs = [trials(indices.correct).reactTimeMS];    % release relative to stimTime on hits
  faTimes = preStimMS(indices.fa);
  faRTs = [trials(indices.fa).reactTimeMS];          % response relative to stimTime on FAs
  missTimes = preStimMS(indices.miss);
  if isempty(hitTimes)
    hitTimes = -10000;
    hitRTs = 0;
  end
  if isempty(faTimes)
     faTimes = -10000;
     faRTs = 0;
  end
  if isempty(missTimes)
     missTimes = -10000;
  end
  releaseTimes = [(hitTimes + hitRTs), (faTimes + faRTs)];

  % histc uses edges in a ridiculous way.  The first and last edge limit the range, but the bin edges used in between
  % aren't the remaining values.  They are the spots halfway between the remaining values.  Absurd.  This is presumably
  % why Matlab recommends histcounts.  But we don't have that function on all our versions of Matlab.  We set the bin
  % edges up so that there are 2 small bins on either side of the range, which we will merge at the end of the counting.
  subplot(dParams.plotLayout{:}, plotIndex);
  cla;
  hold on;
  trialBins = max(10, length(trials) / 20);                    	% 10 or more bins
  set(gca, 'XLim', [1 trialBins]);
  if isfield(file, 'responseLimitMS')
    timeRangeMS = file.preStimMaxMS + file.responseLimitMS;
    minStimXPos = file.preStimMinMS / timeRangeMS * trialBins + 1;
    maxStimXPos = file.preStimMaxMS / timeRangeMS * trialBins + 1;
    if maxStimXPos > minStimXPos
      set(gca,'XTick', [1, minStimXPos, maxStimXPos, trialBins]);
      set(gca,'XTickLabel',{'0', sprintf('%d', file.preStimMinMS), sprintf('%d', file.preStimMaxMS), ...
            sprintf('%d', timeRangeMS)});
    else
      set(gca,'XTick', [1, maxStimXPos, trialBins]);
      set(gca,'XTickLabel',{'0', sprintf('%d', file.preStimMaxMS), sprintf('%d', timeRangeMS)});
    end
  else
    timeRangeMS = max(file.preStimMaxMS, 1);
    minStimXPos = file.preStimMinMS / timeRangeMS * trialBins + 1;
    maxStimXPos = file.preStimMaxMS / timeRangeMS * trialBins + 1;
    set(gca,'XTick', [1, trialBins]);
    set(gca,'XTickLabel',{'0', sprintf('%d', timeRangeMS)});  
  end
  binWidth = timeRangeMS / trialBins;                            % binWidth in MS
  edges = [0.0, linspace(0.5 * binWidth, timeRangeMS - 0.5 * binWidth, trialBins), timeRangeMS];
  counts = hist(hitTimes, edges);
  nHits = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
  counts = hist(faTimes, edges);
  nFAs = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
  counts = hist(missTimes, edges);
  nMisses = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
  nTotal = nHits + nFAs + nMisses;
  counts = hist(releaseTimes, edges);
  nReleases = [counts(1) + counts(2), counts(3:end-2), counts(end-1) + counts(end)];
  if sum(nTotal) < 1                                     % no counts yet?
     return;
  end

  pHits = nHits ./ nTotal;                               % proportions of hits, FAs and misses
  pFAs = nFAs ./ nTotal;
  pMisses = nMisses ./ nTotal;
  pReleases = nReleases ./ (max(nReleases) * 1.25);
  xlabel('Time of stimulus onset (ms)');
  ylabel('Proportion of trials');
  set(gca, 'YLim', [0 1]);
  set(gca,'YTick', [0 1]);
  title('Outcomes Over Trial');

  plot(pHits, 'color', [0.0, 0.7, 0.0], 'lineWidth', 1);
  plot(pFAs, 'color', [0.9, 0.0, 0.0], 'lineWidth', 1);
  plot(pMisses, 'color', [0.6, 0.4, 0.2], 'lineWidth', 1);
  plot(pReleases, 'color', [0.6, 0.6, 0.6], 'lineWidth', 1);
  if maxStimXPos > minStimXPos
    plot([minStimXPos, minStimXPos], [0, 1], 'k:');
  end
  plot([maxStimXPos, maxStimXPos], [0, 1], 'k:');
end
