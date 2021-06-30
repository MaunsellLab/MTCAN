function [] = plotOneHist(i, spikeHist, stimDurMS, histPrePostMS)

  subplot(6, 4, 12+i);
  plot(smooth(spikeHist, min(0.1, 3000 / sum(spikeHist))));
  xticks([histPrePostMS, histPrePostMS + stimDurMS])
  xticklabels({'0', sprintf('%d', stimDurMS)});
  % title([sprintf('%d', (d - 1) * 30), char(176)]); 
end