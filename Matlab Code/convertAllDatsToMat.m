function convertAllDatsToMat()

% Function will convert unconverted .dat files in Lablib Data folder 
% to .mat file and save it in Matlab Data folder

myDir = uigetdir; % sets search dir
myFiles = dir(fullfile(myDir, '*.dat')); 
converted_cnt = 0;

% scans through files in LLDir, and converts
% them if they're not converted
for k = 1:length(myFiles)
  fileNameMat = strrep(myFiles(k).name, '.dat', '.mat');
  cd '../Matlab Data'/
  if ~isfile(fileNameMat)
    convertGRF(myFiles(k).name)
    converted_cnt = converted_cnt + 1;
  end    
end
if converted_cnt == 0
    disp('All files already converted')
end