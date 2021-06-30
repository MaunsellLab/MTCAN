function [file, trials] = convertGRF(fileName, varargin)
% Function loadGRF.m loads Knot (*.dat) files for protocol 'GaborRFMap'into matlab *.mat with
% all the relavant experimental parameters
% Requires readLLFile.m

% Inputs:
%          direct          : file directory
%          datFile         : Name of the *.dat file
%          block           : Block identity with #trials, a matrix of 2
%                            columns. 1st column is block identity, 2nd
%                            column if total no. of trials in every block.
%
%                            example
%                            blocks =[...
%                                    %block    trials#
%                                     1        4*(25+10)
%                                     2        4*(25+10)
%                                     3        4*(25+10)
%                                     4        4*(25+10)
%                                    ];
%                            Where #trials in the valid location (attended)   = 25*4
%                                  #trials in the valid location (attended)   = 10*4
%
%
% Outputs:
%          trialAll          : structure containing all trials variables (excluding instruction/ decoy trials)
%
% Example:
%
%  clear all; close all; clc;
%  direct = '/Data/JohnMaunsellLab/Poseidon_SC_LC/Poseidon_SC_2017-07-27/';
%  datFile = 'Poseidon_2017-07-27_SD4_sel01';
%
%  Blk =[...
%         1  4*(25+10)
%         2  4*(25+10)
%         3  4*(25+10)
%         4  4*(25+10)
%       ];
%  trialAll = loadSD4(direct,datFile,'block',Blk);
%
%
%%-------------------------------------------------------------------------------------------
%
% inputs:

% Blk = getVararg(varargin, 'block');


%% load all the data from all the trials
dirs = targetDirectories();
dataName = strcat(dirs.LLDir, fileName);
file = readLLFile('i', dataName);
numTrials = file.numberOfTrials;
% for t = 1:header.numberOfTrials
% % for t = header.numberOfTrials:-1:1
%   theStruct = readLLFile('t', t);
%   trials(t) = theStruct;
% end

% reformat the relavant data, delete some of the irrelevant data
% for t = 1:length(trials)
%   trials(t).instructTrial = trials(t).trial.data.instructTrial;
%   trials(t).catchTrial = trials(t).trial.data.catchTrial;
%   trials(t).numStim = trials(t).trial.data.numStim;
%   trials(t).targetIndex = trials(t).trial.data.targetIndex;
%   trials(t).targetOnTimeMS = trials(t).trial.data.targetOnTimeMS;
%   trials(t).orientationChangeIndex = trials(t).trial.data.orientationChangeIndex;
%   trials(t).orientationChangeDeg = trials(t).trial.data.orientationChangeDeg;
%   spikeChannels = [trials(27).spike.data.channel];
%   spikeTimes = [trials(27).spike.data.time];
%   trials(t).spike0 = spikeTimes(spikeChannels == 0);
%   trials(t).spike1 = spikeTimes(spikeChannels == 1);
% end
%  clear trials(t).trials;
%  clear trials(t).spike

% trialAll = struct;
% trialAll.header = readLLFile('i', file);
% numTrials = trialAll.header.numberOfTrials;

% Check 1: frame rate at 100 Hz?
% if trialAll.header.displayCalibration.data.frameRateHz ~= 100
%   error('frame rate ~= 100 Hz')
% end

% Initialization
% The trial struct we create will be indexed, so it needs to have the same entries in each entry.  Because different
% trials may be missing some of the entries, we need to initialize them all up front, and load them (or not) with
% each trial individually.  The initialization is a bit differeent for data that is numerical, cell array, etc.

% Numerical data

for dataTypes = {'trialStart', 'trialEnd', 'trialCertify', 'eotCode', 'fixOn', 'fixate', 'instructTrial', ...
      'catchTrial', 'numStim', 'numMap0Stim', 'numMap1Stim', 'photodiodeTime', 'targetIndex', 'targetOnTimeMS', ...
      'oriChangeIndex', 'oriChangeDeg', 'time', 'saccade'}
  dataT = dataTypes{:};
  trials.(dataT) = nan(numTrials, 1);
end

% 1-B: cell array data
%stimulus (sample/test) orientation, on time, and off time
  trials.spike0 = cell(numTrials, 1);
  trials.spike1 = cell(numTrials, 2);
% trialAll.stimOri = cell(numTrials,1);
% [trialAll.stimOri{:}] = deal(nan(2,2));
% trialAll.pulseTrainData = cell(numTrials,1);
% [trialAll.pulseTrainData{:}] = deal(nan(1,1));

% 1-C: eye positions and pupil size
for sides = {'LEye', 'REye'}
  side = sides{:};
  % calibration
  trials.eyeCal.(side).m = cell(numTrials, 1);
  trials.eyeCal.(side).t = cell(numTrials, 1);
  [trials.eyeCal.(side).m{:}] = deal(nan(2, 2));
  [trials.eyeCal.(side).t{:}] = deal(nan(1, 2));
  % raw positions
  trials.eyeRaw.(side) = cell(numTrials, 1);
  [trials.eyeRaw.(side){:}] = deal(nan(1, 2));
  % pupil size
  trials.eyeSize.(side) = cell(numTrials, 1);
  [trials.eyeSize.(side){:}] = deal(nan(1, 2));
  % calibrated positions
%   for periods = {'fixation', 'sample', 'interstim', 'test', 'targetSpotDelay','targetSpot'}
%     per = periods{:};
%     trials.eyePos.(side).(per) = cell(numTrials,1);
%   end
end

hWait = waitbar(0, '', 'name', sprintf('Converting %s...', fileName));

% 2: trial data-------------
  for tIndex = numTrials:-1:1
    waitbar((numTrials - tIndex) / numTrials, hWait, sprintf('Loading trial %d of %d', numTrials - tIndex, numTrials));
    trial = readLLFile('t', tIndex);
    %2-A: time stamps
    trials(tIndex).trialStart = trial.trialStart.timeMS;            % trial start timeStamp
    trials(tIndex).fixOn = trial.fixOn.timeMS;                      % fixation point on timeStamp
    if isfield(trial, 'fixate')
      trials(tIndex).fixate = trial.fixate.timeMS;                  % fixate start
    end
    if isfield(trial, 'saccade')
      trials(tIndex).saccade = trial.saccade.timeMS;                % saccade time stamp
    end
    trials(tIndex).trialEnd = trial.trialEnd.timeMS;                 % trial end timeStamp
    trials(tIndex).eotCode = trial.trialEnd.data;
    if isfield(trial, 'numMap0Stim')
      trials(tIndex).numMap0Stim = trial.numMap0Stim.data;
    else
      trials(tIndex).numMap0Stim = 0;
    end
    if isfield(trial, 'numMap1Stim')
      trials(tIndex).numMap1Stim = trial.numMap1Stim.data;
    else
      trials(tIndex).numMap1Stim = 0;
    end
    if isfield(trial, 'trial')
      trials(tIndex).instructTrial = trial.trial.data.instructTrial;
      trials(tIndex).catchTrial = trial.trial.data.catchTrial;
      trials(tIndex).numTaskStim = trial.trial.data.numTaskStim;
      trials(tIndex).targetIndex = trial.trial.data.targetIndex;
      trials(tIndex).targetOnTimeMS = trial.trial.data.targetOnTimeMS;
      trials(tIndex).orientationChangeIndex = trial.trial.data.orientationChangeIndex;
      trials(tIndex).orientationChangeDeg = trial.trial.data.orientationChangeDeg;
    end
    if isfield(trial, 'photodiodeTime')
      trials(tIndex).photodiodeTime = trial.photodiodeTime.data;
    end
    if isfield(trial, 'spike')
      spikeChannels = [trial.spike.data.channel];
      spikeTimes = [trial.spike.data.time];
      trials(tIndex).spike0 = spikeTimes(spikeChannels == 0);
      trials(tIndex).spike1 = spikeTimes(spikeChannels == 1);
    end
    if isfield(trial, 'stimDesc')
      gaborIndices  = [trial.stimDesc.data.gaborIndex];
      trials(tIndex).taskStimDesc = trial.stimDesc.data(gaborIndices == 0);
      trials(tIndex).map0StimDesc = trial.stimDesc.data(gaborIndices == 1);
      trials(tIndex).map1StimDesc = trial.stimDesc.data(gaborIndices == 2);
    end
    %2-C: stimulus parameter: orientation
%     stimOri = struct2cell(trial.stimDesc.data);                        % stimulus orientations
%     stimOri = squeeze(stimOri(4,:,:));                                 %4th field is 'directionDeg'
%     stimOri = cell2mat(cellfun(@transpose, stimOri, 'uni', 0));
%     trials.stimOri{tIndex} = stimOri;                              % [sample [L R]
%     %  test  ]
%     trials.oriDiff(tIndex,:) = [trial.trial.data.trialSide(1).oriChangeDeg trial.trial.data.trialSide(2).oriChangeDeg];
    
    %2-D: trial type
%     if isfield(trial, 'blockStatus')
%       trials.block(tIndex) = trial.blockStatus.data.blocksDone; % [0, 1, 2,...]
%     end
%     if isfield(trial.trial.data, 'decoyTrial')                         % instruct (to take into account of decoyTrial)
%       trials.instruct(tIndex) = trial.trial.data.instructTrial + trial.trial.data.decoyTrial;
%     else
%       trials.instruct(tIndex) = trial.trial.data.instructTrial;
%     end
%     trials.attendLoc(tIndex) = trial.trial.data.attendLoc;    % [0 (left), 1 (right)]
%     trials.testLoc(tIndex) = trial.trial.data.correctLoc;     % location of test stimulus
%     trials.validTrial(tIndex) = trial.trial.data.validTrial;  % whether test location is on the valid side [0, 1]
%     trials.match(tIndex) = trial.trial.data.nonTargetTrial;   % match trial
%     if isfield(trial.trial.data, 'stimTrial')
%       trials.stimTrial(tIndex) = trial.trial.data.stimTrial;    % stimulus trial
%     end
%     if isfield(trial.trial.data, 'laserPowerMW')
%       trials.laserPowerMW(tIndex) = trial.trial.data.laserPowerMW;   % laserPowerMW
%     end
%     
%     if isfield(trial, 'pulseTrainData')
%       trials.pulseTrainData{tIndex} = trial.pulseTrainData.data;   % pulseTrainData
%     end
%     
%     %2-E: response outcome
%     trials.outcome(tIndex) = trial.signalDetectResult.data; % >0 for completed trials
%     trials.eotcode(tIndex) = trial.trialEnd.data;
    
    %2-F: Eye position
%     if isfield(trial, 'eyeLXData')
%       % calibration
%       trials.eyeCal.LEye.m{tIndex} = [trial.eyeLeftCalibrationData.data.cal.m11 trial.eyeLeftCalibrationData.data.cal.m12;...
%         trial.eyeLeftCalibrationData.data.cal.m21 trial.eyeLeftCalibrationData.data.cal.m22]; %left eye
%       trials.eyeCal.LEye.t{tIndex} = [trial.eyeLeftCalibrationData.data.cal.tX trial.eyeLeftCalibrationData.data.cal.tY];
%       
%       trials.eyeCal.REye.m{tIndex} = [trial.eyeRightCalibrationData.data.cal.m11 trial.eyeRightCalibrationData.data.cal.m12;...
%         trial.eyeRightCalibrationData.data.cal.m21 trial.eyeRightCalibrationData.data.cal.m22] ;%Right eye
%       trials.eyeCal.REye.t{tIndex} = [trial.eyeRightCalibrationData.data.cal.tX trial.eyeRightCalibrationData.data.cal.tY];
%       
%       %raw eye position
%       minLengLeft = min([length(trial.eyeLXData.data), length(trial.eyeLYData.data)]);    %min length of eye data
%       minLengRight = min([length(trial.eyeRXData.data), length(trial.eyeRYData.data)]);
%       
%       trials.eyeRaw.LEye{tIndex} = [trial.eyeLXData.data(1:minLengLeft),trial.eyeLYData.data(1:minLengLeft)]; %left eye
%       trials.eyeRaw.REye{tIndex} = [trial.eyeRXData.data(1:minLengRight),trial.eyeRYData.data(1:minLengRight)]; %left eye
%       trials.eyeSize.LEye{tIndex} = trial.eyeLPData.data; %left eye
%       trials.eyeSize.REye{tIndex} = trial.eyeRPData.data; %right eye
%       
%       %timeStamps of eyeLink data recording start
%       trials.eyeStart(tIndex) = trial.eyeLXData.timeMS(1); %left eye
%     end %Eye position
%   end
  % trial.signalDetectResult.data
  %clear trial
  end % trialInd
  delete(hWait)
  [~, nameOnly, ~] = fileparts(fileName);
  save(strcat(dirs.matlabDir, nameOnly, '.mat'), 'file', 'trials');

  %3-A calibrate eye positions
%   for sides = {'LEye', 'REye'}
%     side = sides{:};
%     trials.eyePos.(side).all = cellfun( @(eyeRaw, m, t) ...
%       (eyeRaw * m + ones(size(eyeRaw,1),1)*t), ...
%       trials.eyeRaw.(side), ...
%       trials.eyeCal.(side).m, ...
%       trials.eyeCal.(side).t, 'UniformOutput', false);
%   end% sides

  %3-B: Updating block status
%   indBlk = unique(Blk);
%   nMissing =[];
%   BlkT =[];
%   for ib=1:length(indBlk)
%     indB = find(Blk==indBlk(ib));
%     BlkT{ib} = Blk(indB); %#ok<AGROW>
%     nMissing(ib) = length(find(diff(trials.attendLoc(indB))~=0)); %#ok<AGROW>
%     BlkT{ib} = BlkT{ib}(1:end-nMissing(ib)); %#ok<AGROW>
%     Blk=[cell2mat(BlkT') ;Blk(indB(end)+1:end)];
%   end

% trialAll.block = NaN(length(trialAll.block),1) ;
% trialAll.block(1:length(Blk)) = Blk;


%% optional inputs
%   function a = getVararg(v, param)
%     ind = find(strcmpi(v, param)) + 1;
%     if isempty(ind)
%       a = [];
%     else
%       a = v{ind};
%     end
%   end % getVararg

end % loadGRF.m