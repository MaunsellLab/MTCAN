function outcomesOverDay(dParams, plotIndex, trials, indices)

% Create plot showing the fractions of hits, misses and FAs over all the trials in the session

numTrials = length(trials);
hits = zeros(1, numTrials);
fas = zeros(1, numTrials);
misses = zeros(1, numTrials);
hits(indices.correct) = 1;
fas(indices.fa) = 1;
misses(indices.miss) = 1;

h = subplot(dParams.plotLayout{:}, plotIndex);
cla reset;
yyaxis left
plot(smooth(hits, ceil(numTrials / 5), 'lowess'), 'color', [0.0, 0.7, 0.0], 'lineWidth', 1);
hold on;
plot(smooth(fas, ceil(numTrials / 5), 'lowess'), 'color', [0.9, 0.0, 0.0], 'lineStyle', '-', 'lineWidth', 1);
plot(smooth(misses, ceil(numTrials / 5), 'lowess'), 'color', [0.6, 0.4, 0.2], 'lineStyle', '-', 'lineWidth', 1);
xlabel('Trials');
ylabel('Proportion of trials');
set(gca, 'YColor', 'k');
set(gca, 'YLim', [0 1]);
set(gca,'YTick', [0, 1]);
title('Outcomes Over Day');
hold off;

yyaxis right
if isfield(trials, 'acquireTimeS')
  plot(smooth([trials(:).acquireTimeS], ceil(numTrials / 5), 'lowess'), ...
       'color', [0.6, 0.6, 0.6], 'lineStyle', '-', 'lineWidth', 1);
end
set(h, 'YColor', 'k');
ylabel('Settle time (s)');
a = axis;
axis([0, numTrials, 0, a(4)])
yh = get(gca, 'ylabel');     % handle to the label object
p = get(yh, 'position');     % get the current position property
p(1) = p(1) * 0.95;
set(yh, 'position', p, 'rotation', -90, 'verticalAlignment', 'middle');      % set the new position, rotate label
set(gca,'YTick', [0, a(4)]);

end
