
% create a list of the first trial sequence
% go through second trial sequence - check to see if sequence is same
%   - if sequence is not same then return error message
%   - if second sequence is greater length, add extra stims to first
%   sequence list
% go through rest of the trials
% if trials not the sequence then return error message

% % for i = 0:2
% %     for j = 0:2
% %         my_lst(end+1,:) = [i,j]; 
% %     end
% % end

lstLoc0 = cell(1,0);
lstLoc1 = cell(1,0);
lstLoc2 = cell(1,0);
lstLoc3 = cell(1,0);

for t = 1:length(trials)
    cur_trial = trials(t);
    if cur_trial.trial.data.instructTrial == 0
        
        % checks for the stimulus in which the target appears
        for d = length(cur_trial.stimDesc.data)
            if ismember(2, cur_trial.stimDesc.data(d).listTypes) == 1
            targetOnsetStim = d;
            end
        end
        
        % checks for attended location and operates on each location 
        % independently 
        if cur_trial.trial.data.attendLoc == 0
            for d = 1:targetOnsetStim-1
                stimSeq = reshape(cur_trial.stimDesc.data(d).stimTypes, [1,4]);
                if length(lstLoc0) < 9
                    lstLoc0{:,end+1} = stimSeq(3:4);
                    currenStim = stimSeq(3:4)
                else
                    %check to see if stim_seq follows same order
                   
                end
            end
        elseif cur_trial.trial.data.attendLoc == 1
            % do another operation
        elseif cur_trial.trial.data.attendLoc == 2
            % do another operation
        else
            % final operation
        end
        
        % check stim sequence
        for d = 1:targetOnsetstim-1
            stimSeq = reshape(cur_trial.stimDesc.data(d).stimTypes, [1,4]);
        end    
    end
end
