function trialAll = loadSD4(direct, datFile, varargin)
% Functiom loadSD4.m loads raw *.dat (knot) files for protocol 'SignalDetection4'into matlab *.mat with
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
         
      Blk = getVararg(varargin, 'block');

  
%% loading trials----------------------------------------------------------------------------
      file=strcat(direct,datFile,'.dat');
      hw = waitbar(0,('Loading trials...'),'position',[20 600 370 60]);
      set(hw,'Name','Please wait');
   
      trialAll = struct;
      trialAll.header = readLLFile('i', file);
      numTrials = trialAll.header.numberOfTrials;
   
      % Check 1: frame rate at 100 Hz?
     if trialAll.header.displayCalibration.data.frameRateHz ~= 100
        error('frame rate ~= 100 Hz')
     end
 
     %Initialization-------------
     % 1-A: numeric data
     for dataTypes = {'block','instruct', 'attendLoc','validTrial','match','testLoc','eotcode','outcome','reward','deliveredReward','stimTrial','laserPowerMW','pulseTrainData',...
                     'trialStart', 'eyeStart','fixOn','fixate','sampleOn','sampleOff','testOn','testOff','targetOn','fixOff','react','saccade','juiceOn','trialEnd',...
                     'fixInterval','sampleInterval','interstim','testInterval', 'targetSpotDelay','fixOffDelay','tooFast'}
         dataT = dataTypes{:};
         trialAll.(dataT) = nan(numTrials, 1);
     end
     trialAll.oriDiff=nan(numTrials,2);
     
     % 1-B: cell array data
     %stimulus (sample/test) orientation, on time, and off time
     trialAll.stimOri = cell(numTrials,1);
     [trialAll.stimOri{:}] = deal(nan(2,2));
     trialAll.pulseTrainData = cell(numTrials,1);
     [trialAll.pulseTrainData{:}] = deal(nan(1,1));
     
     % 1-C: eye positions and pupil size
     for sides = {'LEye', 'REye'}; side = sides{:};
         % calibration
         trialAll.eyeCal.(side).m = cell(numTrials,1);
         trialAll.eyeCal.(side).t = cell(numTrials,1);
         [trialAll.eyeCal.(side).m{:}] = deal(nan(2,2));
         [trialAll.eyeCal.(side).t{:}] = deal(nan(1,2));
         % raw positions
         trialAll.eyeRaw.(side) = cell(numTrials,1);
         [trialAll.eyeRaw.(side){:}] = deal(nan(1,2));
         % pupil size
         trialAll.eyeSize.(side) = cell(numTrials,1);
         [trialAll.eyeSize.(side){:}] = deal(nan(1,2));
         % calibrated positions
         for periods = {'fixation', 'sample', 'interstim', 'test', 'targetSpotDelay','targetSpot'} 
             per = periods{:};       
             trialAll.eyePos.(side).(per) = cell(numTrials,1);
         end
     end
     
     % 2: trial data-------------
     trialInd =0;
     for j = 1:numTrials
         trial = readLLFile('t',j);
         if trial.signalDetectResult.data >= 0                                  %Required trial types (completed trials)
            trialInd = trialInd+1;
            trialAll.indices(trialInd) = j;
            
            %2-A: time stamps
            trialAll.trialStart(trialInd) = trial.trialStart.timeMS;            % trial start timeStamp
            trialAll.fixOn(trialInd) = trial.fixOn.timeMS;                      % fixation point on timeStamp 
            if isfield(trial, 'fixate')
               trialAll.fixate(trialInd) = trial.fixate.timeMS;                 % fixate start
            end
            trialAll.sampleOn(trialInd) = trial.stimulusOn.timeMS(1);           % sample ON timeStamp
            trialAll.testOn(trialInd) = trial.stimulusOn.timeMS(2);             % test ON timeStamp
            trialAll.sampleOff(trialInd) = trial.stimulusOff.timeMS(1);         % sample OFF timeStamp
            trialAll.testOff(trialInd) = trial.stimulusOff.timeMS(2);           % test OFF timeStamp
            if isfield(trial, 'targetSpots')
               trialAll.targetOn(trialInd) = trial.targetSpots.timeMS;           % target on time stamp
            end
            trialAll.fixOff(trialInd) = trial.fixOff.timeMS;                     % fixOff time stamp
            if isfield(trial, 'react')
               trialAll.react(trialInd) = trial.react.timeMS;                    % react time stamp
            end
            if isfield(trial, 'saccade')
               trialAll.saccade(trialInd) = trial.saccade.timeMS;                % saccade time stamp
            end
            if isfield(trial, 'juiceOn')
               trialAll.juiceOn(trialInd) = trial.juiceOn.timeMS(1);                % juiceOn time stamp
            end
            trialAll.trialEnd(trialInd) = trial.trialEnd.timeMS;                 % trial end timeStamp
            
            %2-B: intervals
            trialAll.reward(trialInd)= trial.trial.data.rewardMS;                   %expected reward interval (ms)
            if isfield(trial, 'deliveredRewardMS')
               trialAll.deliveredReward(trialInd)= trial.deliveredRewardMS.data(1); %obtained reward interval (ms)
            else
               trialAll.deliveredReward(trialInd) =0;
            end
            trialAll.tooFast(trialInd) = trial.tooFast.timeMS;                    % too Fast interval (ms)
            trialAll.fixInterval = trialAll.sampleOn - trialAll.fixate;           %fixate interval (ms)
            trialAll.sampleInterval = trialAll.sampleOff - trialAll.sampleOn;     %sample interval (ms)
            trialAll.interstim = trialAll.testOn - trialAll.sampleOff;            %interstim interval (ms)
            trialAll.testInterval = trialAll.testOff - trialAll.testOn;           %test interval (ms)
            trialAll.targetSpotDelay = trialAll.targetOn - trialAll.testOff;      %targetSpot Delay interval (ms)
            trialAll.fixOffDelay = trialAll.fixOff - trialAll.targetOn;           %fixOff Delay interval (ms)
  
            %2-C: stimulus parameter: orientation
            stimOri = struct2cell(trial.stimDesc.data);                        % stimulus orientations
            stimOri = squeeze(stimOri(4,:,:));                                 %4th field is 'directionDeg'
            stimOri = cell2mat(cellfun(@transpose, stimOri, 'uni', 0));
            trialAll.stimOri{trialInd} = stimOri;                              % [sample [L R]
                                                                               %  test  ] 
            trialAll.oriDiff(trialInd,:) = [trial.trial.data.trialSide(1).oriChangeDeg trial.trial.data.trialSide(2).oriChangeDeg];
            
            %2-D: trial type
            if isfield(trial, 'blockStatus')
               trialAll.block(trialInd) = trial.blockStatus.data.blocksDone; % [0, 1, 2,...]
            end
            if isfield(trial.trial.data, 'decoyTrial')                         % instruct (to take into account of decoyTrial)
               trialAll.instruct(trialInd) = trial.trial.data.instructTrial + trial.trial.data.decoyTrial;
            else
               trialAll.instruct(trialInd) = trial.trial.data.instructTrial;
            end
            trialAll.attendLoc(trialInd) = trial.trial.data.attendLoc;    % [0 (left), 1 (right)]
            trialAll.testLoc(trialInd) = trial.trial.data.correctLoc;     % location of test stimulus
            trialAll.validTrial(trialInd) = trial.trial.data.validTrial;  % whether test location is on the valid side [0, 1]
            trialAll.match(trialInd) = trial.trial.data.nonTargetTrial;   % match trial
            if isfield(trial.trial.data, 'stimTrial') 
               trialAll.stimTrial(trialInd) = trial.trial.data.stimTrial;    % stimulus trial
            end
            if isfield(trial.trial.data, 'laserPowerMW')
               trialAll.laserPowerMW(trialInd) = trial.trial.data.laserPowerMW;   % laserPowerMW
            end
            
            if isfield(trial, 'pulseTrainData')
               trialAll.pulseTrainData{trialInd} = trial.pulseTrainData.data;   % pulseTrainData
            end
            
            %2-E: response outcoopme
            trialAll.outcome(trialInd) = trial.signalDetectResult.data; % >0 for completed trials
            trialAll.eotcode(trialInd) = trial.trialEnd.data;
            
            %2-F: Eye position
            if isfield(trial, 'eyeLXData')
               % calibration
               trialAll.eyeCal.LEye.m{trialInd} = [trial.eyeLeftCalibrationData.data.cal.m11 trial.eyeLeftCalibrationData.data.cal.m12;...
                                              trial.eyeLeftCalibrationData.data.cal.m21 trial.eyeLeftCalibrationData.data.cal.m22]; %left eye
               trialAll.eyeCal.LEye.t{trialInd} = [trial.eyeLeftCalibrationData.data.cal.tX trial.eyeLeftCalibrationData.data.cal.tY];

               trialAll.eyeCal.REye.m{trialInd} = [trial.eyeRightCalibrationData.data.cal.m11 trial.eyeRightCalibrationData.data.cal.m12;...
                                              trial.eyeRightCalibrationData.data.cal.m21 trial.eyeRightCalibrationData.data.cal.m22] ;%Right eye
               trialAll.eyeCal.REye.t{trialInd} = [trial.eyeRightCalibrationData.data.cal.tX trial.eyeRightCalibrationData.data.cal.tY];
             
               %raw eye position
               minLengLeft = min([length(trial.eyeLXData.data), length(trial.eyeLYData.data)]);    %min length of eye data
               minLengRight = min([length(trial.eyeRXData.data), length(trial.eyeRYData.data)]);
               
               trialAll.eyeRaw.LEye{trialInd} = [trial.eyeLXData.data(1:minLengLeft),trial.eyeLYData.data(1:minLengLeft)]; %left eye
               trialAll.eyeRaw.REye{trialInd} = [trial.eyeRXData.data(1:minLengRight),trial.eyeRYData.data(1:minLengRight)]; %left eye
               trialAll.eyeSize.LEye{trialInd} = trial.eyeLPData.data; %left eye
               trialAll.eyeSize.REye{trialInd} = trial.eyeRPData.data; %right eye
            
               %timeStamps of eyeLink data recording start
               trialAll.eyeStart(trialInd) = trial.eyeLXData.timeMS(1); %left eye
            end %Eye position
         end% trial.signalDetectResult.data
         waitbar(j/numTrials,hw,strcat('Loading trials...',num2str(j)));
         %clear trial
     end % trialInd
     delete(hw)
   
     %3-A calibrate eye positions
     for sides = {'LEye', 'REye'} 
         side = sides{:};
         trialAll.eyePos.(side).all = cellfun( @(eyeRaw, m, t) ...
                                               (eyeRaw*m + ones(size(eyeRaw,1),1)*t), ...
                                               trialAll.eyeRaw.(side), ...
                                               trialAll.eyeCal.(side).m, ...
                                               trialAll.eyeCal.(side).t, 'UniformOutput', false);  
     end% sides

     %3-B: Updating block status
     indBlk = unique(Blk);
     nMissing =[];
     BlkT =[];
     for ib=1:length(indBlk)
         indB = find(Blk==indBlk(ib));
         BlkT{ib} = Blk(indB);
         nMissing(ib) = length(find(diff(trialAll.attendLoc(indB))~=0));
         BlkT{ib} = BlkT{ib}(1:end-nMissing(ib));
         Blk=[cell2mat(BlkT') ;Blk(indB(end)+1:end)];
     end
     
     trialAll.block = NaN(length(trialAll.block),1) ;
     trialAll.block(1:length(Blk)) = Blk;
     
 
 %% optional inputs
  function a = getVararg(v, param)
        ind = find(strcmpi(v, param)) + 1;
        if isempty(ind)
            a = [];
        else
            a = v{ind};
        end
  end % getVararg

end % loadSD4.m