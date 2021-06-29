function testData()

  fileName = 'testc00.dat';
%   convertGRF(fileName);
  
  dirs = targetDirectories();
  [~, nameOnly, ~] = fileparts(fileName);
  load(strcat(dirs.matlabDir, nameOnly, '.mat'), 'file', 'trials');
  
  if file.mapSettings.data.directionDeg.n > 1
    for c = 0:1
      plotOneChanneldirectionDeg(c, file, trials);
    end
  else
    fprintf('no direction tuning\n');
  end
  
  if file.mapSettings.data.temporalFreqHz.n > 1
    for c = 0:1
      plotOneChanneltemporalFreqHz(c, file, trials);
    end
  else
    fprintf('no speed tuning\n');
  end
  
  if file.mapSettings.data.azimuthDeg.n > 1
    for c = 0:1
      plotOneChannelRFHeatMap(c, file, trials);
    end
  else
    fprintf('no RF spatial map\n');
  end

%   if file.mapSettings.data.directionDeg.n <= 1
%     fprintf('no direction tuning\n');
%     return;
%   end
%   for c = 0:1
%     plotOneChanneldirectionDeg(c, file, trials);
%   end
%   if file.mapSettings.data.temporalFreqHz.n <= 1
%     fprintf('no speed tuning\n');
%     return;
%   end
%   for c = 0:1
%     plotOneChanneltemporalFreqHz(c, file, trials);
%   end
end

