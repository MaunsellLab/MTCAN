%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Text
function drawText(taskName, dParams, file, trials, indices)
  axisHandle = subplot(dParams.plotLayout{:}, 1);						% default axes are 0 to 1
  set(axisHandle, 'Visible', 'off');
  set(axisHandle, 'OuterPosition', [0.02 0.75, 0.25, 0.2]);
  text(0.00, 1.25, taskName, 'FontWeight', 'bold', 'FontSize', 16);
  if (nargin < 3)                                          % display task title only
    return;
  end
  muStr = '\mul';
  if file.subjectNumber == 0 || file.subjectNumber == 999
      text(0.00, 1.11, sprintf('Subject: %d', file.subjectNumber), 'FontSize', 14, 'Color', 'r');
  else
      text(0.00, 1.11, sprintf('Subject: %d', file.subjectNumber), 'FontSize', 14);
  end
  headerText = cell(1, 1);
  if isfield(file, 'startTimeVec')
      headerText{1} = datestr(file.startTimeVec, 'mmmm dd, yyyy HH:MM');
  else
      headerText{1} = '(missing date field)';
  end
  if isfield(trials(:), 'rewardUL')
      headerText{length(headerText) + 1} = sprintf('%.0f minutes, %.0f %s', ...
          (trials(end).totalRunTimeS - trials(1).totalRunTimeS) / 60.0, sum([trials(:).rewardUL]), muStr);
  else
      headerText{length(headerText) + 1} = sprintf('%.0f minutes', ...
          (trials(end).totalRunTimeS - trials(1).totalRunTimeS) / 60.0);
  end
  text(0.00, 1.00, headerText, 'units', 'normalized', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');

  % Display trial outcomes in a table
  numCorrect = sum(indices.correct);
  numFailed = sum(indices.miss);
  numEarly = sum(indices.fa);
  totalTrials = numCorrect + numFailed + numEarly;
  t2H(1) = text(0.00, 0.77, {'Trials:', 'Correct:', 'Failed:', 'Early:'});
  t2H(2) = text(0.35, 0.77, {sprintf('%d', file.trials), sprintf('%d', numCorrect), sprintf('%d', numFailed), ...
          sprintf('%d', numEarly)});
  t2H(3) = text(0.60, 0.77, {' ', sprintf('%.0f%%', numCorrect / totalTrials * 100.0), ...
          sprintf('%.0f%%', numFailed / totalTrials * 100.0), sprintf('%.0f%%', numEarly / totalTrials * 100.0)});
  set(t2H, 'units', 'normalized', 'verticalAlignment', 'top', 'horizontalAlignment', 'left');
  set(gcf, 'Visible', 'on');
  if (isfield(trials, 'syntheticData') && trials(end).trial.syntheticData == 1) || ...
      (isfield(file, 'syntheticData') &&  file.syntheticData == 1)
    text(0.85, 1.33, 'Synthetic Data', 'FontSize', 16, 'Color', 'r', 'fontWeight', 'bold', 'VerticalAlignment', 'top');
  end
end