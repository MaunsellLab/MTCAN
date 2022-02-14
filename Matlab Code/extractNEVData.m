function [taskNames, taskEvents, taskSpikes] = extractNEVData(NEV)

% function extractNEVData:
% Extracts Lablib events and spikes from NEV array returned by readNEV
% Breaks the data into chunks associated with different Lablib task plugins
% Those chunks are returned in taskNames, taskEvents, and taskSpikes, which are the same length,
%  and have one entry each for each run of a Lablib task.
%%-------------------------------------------------------------------------------------------------
%Inputs:    NEV, NEV structure from readNEV
%Outputs:   taskName: cell array of the letter code for each plugin run
%           taskEvents: cell array of cell arrays of event from each plugin run
%           taskSpikes: cell arrray of array of spike events associated with each plugin run
%%-------------------------------------------------------------------------------------------------
addpath('/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NPMK/');
addpath('/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NEV/');
  % for testing
  if nargin < 1
    nevFile = 'testing_220214001';                          % .nev file (no file extension)
    directory = '/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NEV/';  % directory for .nev file
    NEV = readNEV([directory, nevFile, '.nev']);          	                      % read .nev file
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
    index = startIndices(t);                                % event for trial start
    eventCounts = zeros(1, length(codes));
    events(t).trialStart.timeS = timeS(index);
    events(t).trialStart.data = data(index);
    while true
      index = index + 1;                                  % advance to next event
      if index > size(data, 1)                            % end of data, we're done;
        break;
      end
      if codes(index) == endCode                          % end of trial, process & quit            
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
      if ~isfield(events(t), eventName)                   % first occurrence this trial
        events(t).(eventName{1}) = event;
      else                                                % multiple occurrences this trial
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
  taskEnds = [find(diff(taskCodes) ~= 0), length(taskCodes)]; % indices for trial of each task
  taskNames = cell(1, length(taskEnds));
  taskEvents = cell(1, length(taskEnds));
  taskSpikes = cell(1, length(taskEnds));
  taskStart = 1;                                              % start of first task run
  for t = 1:length(taskEnds)
    taskNames{t} = dict.taskNames{taskCodes(taskStart)};   	% name of task
    taskEvents{t} = events(taskStart:taskEnds(t));           % extract the events for this task run
    startTime = taskCodesStruct(taskStart).time;             % startTime of taskRun
    if t < length(taskEnds)
      endTime = taskCodesStruct(taskEnds(t) + 1).time;
    else
      endTime = inf;
    end
    taskSpikes{t} = spikes(spikes(:,3) >= startTime & spikes(:,3) < endTime, :);
    taskStart = taskEnds(t) + 1;                             % advance to start of next task run
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
%     letterCodes = sort(unique(codes));                      % all codes used
%     factors = zeros(1, length(letterCodes));
%     for i = 1:length(letterCodes)
%         code = dec2bin(letterCodes(i));
%         letterCode = [char(bin2dec(code(2:8))), char(bin2dec(code(10:16)))];
% %         index = find(ismember(taskLetters, letterCode), 1);         % get the index for the code
% %         index = find(sum(ismember(taskLetters, letterCode), 2) == 2);
%         for index = 1:size(dict.letters, 1)
%             if dict.letters(index,:) == letterCode
%                 break;
%             end
%         end
%         if index > size(dict.letters, 1)                  	% no match?
%             error('defineEvents: failed to find name for code ''%s''', letterCode);
%         end
%         names(i) = dict.names(index);                      	%#ok<AGROW> % load the name into the dictionary
%         factors(i) = dict.multipliers(index);               % save the sorted factors (temporarily)
%     end
    % convert all data words by bit shifting and multiplying as needed
    for i = 1:length(dict.multipliers)
        if dict.multipliers(i) ~= 1.0
            data(codes == dict.codes(i)) = data(codes == dict.codes(i)) .* dict.multipliers(i);
%             factor = factors(find(letterCodes == codes(i), 1));
%             dataWord = dec2bin(data(i), 16);
%             data(i) = bin2dec(dataWord(2:15)) * factor;
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
    %showTrouble(index, codes, data);                            % for debugging
    if data(index - 1) > 2^15                                   % spurious code word
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
        followers = codes(indices + 1);                         % the set of codes follwing the preceding code
        firstProb = length(find(followers == firstCode)) / length(followers);
        secondProb = length(find(followers == secondCode)) / length(followers);
        if firstProb >= secondProb
            dEvents(index * 2 - 2, :) = [];                     % delete the second
        else
            dEvents(index * 2 - 3, :) = [];                     % delete the first
        end
        fprintf('Warning: Two valid codes in sequence at index %d: ''%s'', ''%s''.', index - 1, firstLetters,...
            secondLetters);
        fprintf(' Corrected by deleting the less probable (%.3f versus %.3f)\n', firstProb, secondProb);
      else
        error('extractNEVData:validateData: Two invalid codes in a row\n');
      end
    else                                                % two data words in a row
      preCode = codes(index - 1);                     % valid code immediately preceding the suspects
      indices = codes(1:index-3) == preCode;    % all occurrences of that code
      followers = data(indices);                      % the set of data following the code
      firstProb = length(find(followers == data(index - 1))) / length(followers);
      secondProb = length(find(followers == codes(index))) / length(followers);
      if firstProb < secondProb
        dEvents(index * 2 - 2, :) = [];             % delete the first data word
      else
        dEvents(index * 2 - 1, :) = [];         	% delete the second data word
      end
      fprintf('Warning: Spurious data word at index %d of %d.', index, length(codes));
      fprintf(' Corrected by deleting the less probable (%.3f versus %.3f)\n', firstProb, secondProb);
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
    error('Error: NEV contains unknown codes.  This is a serious error.');
  end
end

%%
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

%% charsFromCode -- convert an event code into two letters for the code

function [charCode, binaryCode] = convertCode(code)

  binaryCode = dec2bin(code, 16);
  charCode = [char(bin2dec(binaryCode(2:8))), char(bin2dec(binaryCode(10:16)))];
end

