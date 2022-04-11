function [taskNames, taskEvents, taskSpikes] = extractNEVData(NEV)

% function extractNEVData:
% Extracts Lablib events and spikes from NEV array returned by readNEV
% Breaks the data into chunks associated with different Lablib task plugins
%  Distinguishes chucks associated with different task IDs for a given task plugin
% Chunks are returned in taskNames, taskEvents, and taskSpikes, which are the same length,
%  and have one entry each for each run of a Lablib task/ID.
%-------------------------------------------------------------------------------------------------
% Inputs:    NEV, NEV structure from readNEV
% Outputs:   taskName: cell array of the letter code for each plugin run
%           taskEvents: cell array of cell arrays of event from each plugin run
%           taskSpikes: cell arrray of array of spike events associated with each plugin run
%-------------------------------------------------------------------------------------------------
% addpath('/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NPMK/');
% addpath('/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NEV/');
  % for testing
  if nargin < 1
    nevFile = 'testing_220408';                        % .nev file (no file extension)
    directory = '/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NEV/';  % directory for .nev file
    NEV = readNEV([directory, nevFile, '.nev']);              % read .nev file
  end

  % extract events
  dEventIndices = find(NEV(:, 1) == 0);                       % Lablib events are channel 0
  dEvents = NEV(dEventIndices, :);                            % extract Lablib events
  NEV(dEventIndices, :) = [];                                 % remove Lablib events
  spikes = NEV;                                               % all other events from the NEV are spikes

  % validate events, which should alternate between codes (bit 15 set) and data (bit 15 cleared)
  dict = taskDict();                                          % get dictionary for interpreting data
  [data, codes, timeS] = parseEvents(dEvents, dict);          % parse data 
  endCode = dict.codes(ismember(dict.names, 'trialEnd'));     % get the trialEnd code
  startCode = dict.codes(ismember(dict.names, 'trialStart')); % get the trialStart code
  startIndices = find(codes == startCode);                    % get the trial starts
  events = struct([]);
  for t = 1:length(startIndices)                              % for each trial...
    index = startIndices(t);                                  % event for trial start
    eventCounts = zeros(1, length(codes));
    events(t).trialStart.timeS = timeS(index);
    events(t).trialStart.data = data(index);
    while true
      index = index + 1;                                      % advance to next event
      if index > size(data, 1)                                % end of data, we're done;
        break;
      end
      if codes(index) == endCode                              % end of trial, process & go to next trial            
        events(t).trialEnd.timeS = timeS(index);
        events(t).trialEnd.data = data(index);
        break;
      end
      if codes(index) == startCode
        fprintf('Warning: ''trialStart'' with no ''trialEnd''. Trial %d index %d', t, index);
        break;
      end
      eventIndex = find(dict.codes == codes(index), 1);
      eventName = dict.names(eventIndex);
      eventCounts(index) = eventCounts(index) + 1;
      event.time = timeS(index);
      event.data = data(index);
      if ~isfield(events(t), eventName)                       % first occurrence this trial
        events(t).(eventName{1}) = event;
      else                                                    % multiple occurrences this trial
       	events(t).(eventName{1}) = [events(t).(eventName{1}) event];
      end
    end
  end
  for t = length(startIndices):-1:1                           % delete any trials that got no endTrial code
    if isempty(events(t).trialEnd)
      events(t) = [];
    end
  end

  % Now that we've extracted the sequence of events, we need to split it into data from different plugins
  taskCodesStruct = [events(:).taskCode];                   	% all the taskCode structs
  taskCodes = [taskCodesStruct(:).data];                  	  % all the taskCodes
  taskSet = unique(taskCodes);                                % unique codes
  cases = 0;                                                  % number of setting-tasks
  if (isfield(events, 'settingsCode'))
    settingsCodesStruct = [events(:).settingsCode];           % all the settingCode structs
    settingsCodes = [settingsCodesStruct(:).data];            % all the settingCodes
    settingsSet = unique(settingsCodes);                      % unique codes
    taskSettings = cell(length(settingsSet), length(taskSet));  % logical arrays of all possible setting-tasks
    taskNames = cell(1, length(taskSet) * length(settingsSet)); % names of existing setting-tasks, in order
    for t = 1:length(taskSet)
      for s = 1:length(settingsSet)
        taskSettings{s, t} = taskCodes == taskSet(t) & settingsCodes == settingsSet(s);
        if sum(taskSettings{s, t}) > 0                        % existing setting-task, save a name for it
          cases = cases + 1;
          taskNames{cases} = sprintf('%s%d', dict.taskNames{taskSet(t)}, settingsSet(s));
        end
      end
    end
  else
    settingsSet = 0;                                          % no separate settings, just task.
    taskSettings = cell(1, length(taskSet));
    taskNames = cell(1, length(taskSet));
    for t = 1:length(taskSet)
        taskSettings{t} = taskCodes == taskSet(t);
        if sum(taskSettings(t) > 0)                           % existing task, save a name for it
          cases = cases + 1;
          taskNames{cases} = dict.taskNames{taskSet(t)};
        end
    end
  end

  % Now we've sorted out the setting-tasks, bundle the trial events and spikes for each.
  taskEvents = cell(1, cases);
  taskSpikes = cell(1, cases);
  taskIndex = 1;
  for t = 1:length(taskSet)
    for s = 1:length(settingsSet)
      if sum(taskSettings{s, t}) == 0                         % not an existing settings-task combination
        continue;
      end
      taskEvents{taskIndex} = events(taskSettings{s, t});             % extract the events for this settings-task combination
      firstTrialIndex = find(taskSettings{s, t} > 0, 1, 'first');
      lastTrialIndex = find(taskSettings{s, t} > 0, 1, 'last');
      startTime = taskCodesStruct(firstTrialIndex).time;       % startTime of settings-task run
      if lastTrialIndex < length(taskCodesStruct) 
        endTime = taskCodesStruct(lastTrialIndex + 1).time;    % endTime of settings-task run
      else
        endTime = inf;
      end
      taskSpikes{taskIndex} = spikes(spikes(:,3) >= startTime & spikes(:,3) < endTime, :);    % get the settings-task spikes
      taskIndex = taskIndex + 1;
    end
  end
    
%% Saving the NEV file .mat

  clear codes data dict dEventIndices dEvents endCode endTime event eventCounts events eventIndex eventName index NEV ...
    spikes startCode startIndices startTime t taskCodes taskCodesStruct taskEnds taskStart timeS
% sessionId = input('Enter MouseID and Date: ', 's');
% fileName = strcat(sessionId,'_','NEV','.mat');
  fileName = strcat(nevFile,'_','NEV','.mat');
  save([directory fileName]);
end 

%%
function [data, codes, timeS] = parseEvents(dEvents, dict)
% parseEvents: return aligned arrays of defined codes, data and times for the Lablib events

    [codes, data, timeS] = validateData(dEvents, dict);     % get the codes, data and timestamps
    validateCodes(codes, dict);                             % make sure all codes are known
    % convert all data words by bit shifting and multiplying as needed
    for i = 1:length(dict.multipliers)
        if dict.multipliers(i) ~= 1.0
            data(codes == dict.codes(i)) = data(codes == dict.codes(i)) .* dict.multipliers(i);
        end
    end
end

%%
function [codes, data, timeS] = validateData(dEvents, dict)
% validateData: check a stream of Lablib events for validity.  The sequence should be alternating code and data 
% words
  codes = dEvents(1:2:end - 1, 2);
  data = dEvents(2:2:end, 2);
  timeS = dEvents(1:2:end - 1, 3);
  index = find(codes <= 2^15, 1);
  if ~isempty(index)                                              % trouble: code words with MSB set
%     showTrouble(index, codes, data);                             % for debugging
    if (index == 1)                                               % first code is malformed
      codes(1) = [];
      data(1) = []; %#ok<NASGU> 
      timeS(1) = []; %#ok<NASGU> 
      fprintf('Warning: Spurious code word at index 1 of %d\n', 1, length(codes));
      fprintf(' Corrected by deleting code and data''\n');
    elseif index > 1 && data(index - 1) > 2^15                    % spurious code word
      firstCode = codes(index - 1);                               % get the last (presumably) valid code
      [firstLetters, ~] = convertCode(firstCode);
      firstCodeValid = ismember(firstLetters, dict.letters);
      secondCode = data(index - 1);                               % get the following word (presumably a code)
      [secondLetters, ~] = convertCode(secondCode);
      secondCodeValid = ismember(secondLetters, dict.letters);
      if firstCodeValid && ~secondCodeValid 
        fprintf('Warning: Spurious code word at index %d of %d', index, length(codes));
        fprintf(' Corrected by deleting code ''%s''\n', secondLetters);
        dEvents(index * 2 - 2, :) = [];
      elseif ~firstCodeValid && secondCodeValid
        fprintf('Warning: Spurious code word at index %d of %d.', index, length(codes));
        fprintf(' Corrected by deleting code ''%s''\n', firstLetters);
        dEvents(index * 2 - 3, :) = [];
      elseif firstCodeValid && secondCodeValid
        preCode = codes(index - 2);                             % valid code immediately preceding the suspects
        indices = find(codes(1:index-3) == preCode);            % all occurrences of preceding code
        followers = codes(indices + 1);                         % the set of codes following the preceding code
        firstProb = length(find(followers == firstCode)) / length(followers);
        secondProb = length(find(followers == secondCode)) / length(followers);
        if firstProb >= secondProb
            dEvents(index * 2 - 2, :) = [];                     % delete the second
        else
            dEvents(index * 2 - 3, :) = [];                     % delete the first
        end
        if firstCode == secondCode
          fprintf('Warning: Repeated code at index %d: ''%s''.', index - 1, firstLetters);
          fprintf(' Duplicate deleted.\n');
        else
          fprintf('Warning: Two valid codes in sequence at index %d: ''%s'', ''%s''.', index - 1, firstLetters,...
              secondLetters);
          fprintf(' Deleted the less probable (%.3f versus %.3f)\n', firstProb, secondProb);
        end
      else
        error('extractNEVData:validateData: Two invalid codes in a row\n');
      end
    % two data words in a row
    elseif codes(index) == data(index - 1)                      % repeated data word
      dEvents(index * 2 - 2, :) = [];                           % delete the first data word(index - 1)
      fprintf('Warning: Repeated data at index %d: (%d).', index, codes(index));
      fprintf(' Duplicate deleted.\n');
    else
      preCode = codes(index - 1);                               % valid code immediately preceding the suspects
      indices = codes(1:index-3) == preCode;                    % all occurrences of that code
      followers = data(indices);                                % the set of data words following the code
      firstProb = length(find(followers == data(index - 1))) / length(followers);
      secondProb = length(find(followers == codes(index))) / length(followers);
      if firstProb < secondProb
        dEvents(index * 2 - 2, :) = [];                         % delete the first data word
      else
        dEvents(index * 2 - 1, :) = [];         	              % delete the second data word
      end
      fprintf('Warning: Spurious data at index %d.', index);
      fprintf(' Deleted the less probable (%.3f versus %.3f)\n', firstProb, secondProb);
   end
   [codes, data, timeS] = validateData(dEvents, dict);	% reload and check cleaned up data
  end
  index = find(data > 2^15, 1);
  if ~isempty(index)
      error('Data word at %d has bit 15 set.\nThis needs code development to fix it.', index);
  end
end

%% validateCodes
function validateCodes(codes, dict)
  % Make sure that all the codes from the NEV exist in the dictionary
  codesFromNEV = unique(codes);
  knownCodes = ismember(codesFromNEV, dict.codes);
  if sum(knownCodes) ~= length(codesFromNEV)
    unknownCodes = codesFromNEV(~ismember(codesFromNEV, dict.codes));
    fprintf('Error: NEV contains unknown codes: ');
    for c = 1:length(unknownCodes)
      charCode = convertCode(unknownCodes(c));
      fprintf('''%s''', charCode);
    end
    error('\nUpdate taskDict.m?  Quitting.');
  end
end

function showTrouble(troubleIndex, codes, data) %#ok<*DEFNU>
% Debugging information.  Displays the values and bit patterns around a trouble spot in the event sequence.
  offset = 2;
  startIndex = max(1, troubleIndex - offset);
  endIndex = min(length(codes), troubleIndex + offset);
  for index = startIndex:endIndex
    marker = ' ';
    if index == troubleIndex
      marker = '*';
    end
    [charCode, binaryCode] = convertCode(codes(index));
    fprintf('%d:%s %8d ''%s'' %s\n', index, marker, codes(index), charCode, binaryCode);
    [charCode, binaryCode] = convertCode(data(index));
    fprintf('%d:%s %8d ''%s'' %s\n', index, ' ', data(index), charCode, binaryCode);
 end
end

% charsFromCode -- convert an event code into two letters for the code
function [charCode, binaryCode] = convertCode(code)
  binaryCode = dec2bin(code, 16);
  charCode = [char(bin2dec(binaryCode(2:8))), char(bin2dec(binaryCode(10:16)))];
end

