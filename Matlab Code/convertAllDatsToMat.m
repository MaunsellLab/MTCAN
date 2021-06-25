function convertAllDatsToMat()

% Function will convert unconverted .dat files in Lablib Data folder 
% to .mat file and save it in Matlab Data folder

myDir = uigetdir; % sets search dir
myFiles = dir(fullfile(myDir, '*.dat')); % gets all .dat files into struct

for k = 1:length(myFiles)
  fileName = erase(myFiles(k).name, '.dat');
  fileNameMat = append(fileName, '.mat');
  cd ..
  cd 'Matlab Data'/
  if isfile(fileNameMat)
    % do nothing
  else
    cd ..
    cd 'Matlab Code'/
    convertGRF(myFiles(k).name)
  end    
end
end