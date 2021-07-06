function cumulativeTime(dParams, subplotIndex, trials)
%
% Display the work efficiency of the subject as a percent of maximum efficiency
%
  minRunTimeS = [trials(:).minRunTimeS] - trials(1).minRunTimeS(1);
  totalRunTimeS = [trials(:).totalRunTimeS] - trials(1).totalRunTimeS(1);
  axisHandle = subplot(dParams.plotLayout{:}, subplotIndex);
  cla;
  hold off;
  plot(totalRunTimeS / 60.0, minRunTimeS / 60.0, 'b-', 'linewidth', 1.0);
  hold on
  if length(trials) > 10
    text(0.05, 0.95, sprintf('%.0f trials per hour', 3600.0 * length(trials) / totalRunTimeS(end)), ...
       'units', 'normalized', 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
  end
  xlabel('Actual work time (min)');
  ylabel('Ideal work time (min)');
  a = axis;
  limit = max(a(2), a(4));
  set(axisHandle, 'xlim', [0, limit], 'ylim', [0, limit]);
  plot([0, limit], [0, limit], 'k--');
  efficiency =  mean(minRunTimeS) / mean(totalRunTimeS);
  plot([0, limit], [0, limit * efficiency], 'k:');
  title(sprintf('Work Efficiency %.0f%%', efficiency * 100.0));
end
