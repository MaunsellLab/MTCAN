clear

baseFileName = 'testing_220408';

directory = 'C:\Users\Maunsell Lab\Documents\Chery\cheryKSTest3\lablib';

fileNameMTC = [baseFileName '_MTC.dat'];
fileNameMTN = [baseFileName '_MTN.dat'];
fileNameRFGRF = [baseFileName '_RF_GRF.dat'];
fileNameDirGRF = [baseFileName '_Dir_GRF.dat'];
fileNameSpeedGRF = [baseFileName '_Speed_GRF.dat'];

if isfile(fileNameRFGRF)
% to convert GRF to same format as MTC/MTN
     header = readLLFile('i', fileNameRFGRF);
     trials = cell(1,header.numberOfTrials);
     for i = 1:header.numberOfTrials
         trials{i} = readLLFile('t', i);
     end
     
     fileNameRFGRF = strrep(fileNameRFGRF, '.dat', '.mat');
     save(fileNameRFGRF, 'trials', 'header')
end

if isfile(fileNameDirGRF)
% to convert GRF to same format as MTC/MTN
     header = readLLFile('i', fileNameDirGRF);
     trials = cell(1,header.numberOfTrials);
     for i = 1:header.numberOfTrials
         trials{i} = readLLFile('t', i);
     end
     
     fileNameDirGRF = strrep(fileNameDirGRF, '.dat', '.mat');
     save(fileNameDirGRF, 'trials', 'header')
end

if isfile(fileNameSpeedGRF)
% to convert GRF to same format as MTC/MTN
     header = readLLFile('i', fileNameSpeedGRF);
     trials = cell(1,header.numberOfTrials);
     for i = 1:header.numberOfTrials
         trials{i} = readLLFile('t', i);
     end
     
     fileNameSpeedGRF = strrep(fileNameSpeedGRF, '.dat', '.mat');
     save(fileNameSpeedGRF, 'trials', 'header')
end

if isfile(fileNameMTC)
    header = readLLFile('i', fileNameMTC);
    trials = cell(1,header.numberOfTrials);
    for i = 1:header.numberOfTrials
        trials{i} = readLLFile('t', i);
    end
    
    fileNameMTC = strrep(fileNameMTC, '.dat', '.mat');
    save(fileNameMTC, 'trials', 'header')
end

if isfile(fileNameMTN)
    header = readLLFile('i', fileNameMTN);
    trials = cell(1,header.numberOfTrials);
    for i = 1:header.numberOfTrials
        trials{i} = readLLFile('t', i);
    end
    
    fileNameMTN = strrep(fileNameMTN, '.dat', '.mat');
    save(fileNameMTN, 'trials', 'header')
end


% 
% if isfile(fileNameGRF)
% %     to convert GRF to same format as MTC/MTN
% %     header = readLLFile('i', fileNameGRF);
% %     trials = cell(1,header.numberOfTrials);
% %     for i = 1:header.numberOfTrials
% %         trials{i} = readLLFile('t', i);
% %     end
% %     
% %     cd '../Matlab Data'/
% %     fileNameGRF = strrep(fileNameGRF, '.dat', '.mat');
% %     save(fileNameGRF, 'trials', 'header')
%     convertGRF(fileNameGRF);
%     cd '..'/'Lablib Data'/
% end


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

