function directories = targetDirectories()

% rootDir = '/Users/Shared/Data/MTCAN/';
rootDir = '/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/';

% General path accessing for scripts/data 
% Relies on folder organization of: MTCAN/
%  - Matlab Data/
%  - Matlab Code/
%  - Lablib Data/
%  - BlackRock Data/
%  - Figures/
% cur_location = mfilename('fullpath');
% cur_location_steps = split(cur_location, '/');
% rootDir = strcat(fullfile(steps_to_cur_location{2:end-2}), '/'); 


directories.LLDir = [rootDir, 'Lablib Data/'];
directories.BRDir = [rootDir, 'BlackRock Data/'];
directories.matlabDir = [rootDir, 'Matlab Data/'];
directories.figureDir = [rootDir, 'Figures/'];

end
