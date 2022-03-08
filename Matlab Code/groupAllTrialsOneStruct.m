clear

baseFileName = 'testing_220304';

directory = '/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/Lablib Data/';

fileNameMTC = [baseFileName '_MTC.dat'];
fileNameMTN = [baseFileName '_MTN.dat'];
fileNameGRF = [baseFileName '_GRF.dat'];

if isfile(fileNameGRF)
%     to convert GRF to same format as MTC/MTN
%     header = readLLFile('i', fileNameGRF);
%     trials = cell(1,header.numberOfTrials);
%     for i = 1:header.numberOfTrials
%         trials{i} = readLLFile('t', i);
%     end
%     
%     cd '../Matlab Data'/
%     fileNameGRF = strrep(fileNameGRF, '.dat', '.mat');
%     save(fileNameGRF, 'trials', 'header')
    convertGRF(fileNameGRF);
    cd '..'/'Lablib Data'/
end

if isfile(fileNameMTC)
    header = readLLFile('i', fileNameMTC);
    trials = cell(1,header.numberOfTrials);
    for i = 1:header.numberOfTrials
        trials{i} = readLLFile('t', i);
    end
    
    cd '../Matlab Data'/
    fileNameMTC = strrep(fileNameMTC, '.dat', '.mat');
    save(fileNameMTC, 'trials', 'header')
    cd '..'/'Lablib Data'/
end

if isfile(fileNameMTN)
    header = readLLFile('i', fileNameMTN);
    trials = cell(1,header.numberOfTrials);
    for i = 1:header.numberOfTrials
        trials{i} = readLLFile('t', i);
    end
    
    cd '../Matlab Data'/
    fileNameMTN = strrep(fileNameMTN, '.dat', '.mat');
    save(fileNameMTN, 'trials', 'header')
    cd '..'/'Lablib Data'/
end




% fileName = 'testing_220214_MTC.dat';
% 
% header = readLLFile('i', fileName);
% 
% trials = {};
% for i = 1:header.numberOfTrials
%     trials{i} = readLLFile('t', i);
% end
% 
% cd '../Matlab Data'/
% fileName = strrep(fileName, '.dat', '.mat');
% save(fileName, 'trials', 'header')

% use code below to save as matlab 7.3 version
% save(fileName, 'trials', 'header', '-v7.3');


% if contains(fileName, 'MTC')
%     trialsMTC = trials;
%     cd '../Matlab Data'/
%     fileName = strrep(fileName, '.dat', '.mat');
%     save(fileName, 'trialsMTC', 'header')
% elseif contains(fileName, 'MTN')
%     trialsMTN = trials;
%     cd '../Matlab Data'/
%     fileName = strrep(fileName, '.dat', '.mat');
%     save(fileName, 'trialsMTN', 'header')
% elseif contains(fileName, 'GRF')
%     trialsGRF = trials;
%     cd '../Matlab Data'/
%     fileName = strrep(fileName, '.dat', '.mat');
%     save(fileName, 'trialsGRF', 'header')
% end
% 
% clear trials

