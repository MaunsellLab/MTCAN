function [inited, dParams, file] = checkInitialization(dParams, file)
% Check that initialization has been completed for online display of results

  % Check that the figure has been set up.  Set if up if it hasn't been prepared yet
  h = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', 1);
  if isempty(h) || isempty(dParams) || ~strcmp(get(h, 'type'), 'figure') || nargin == 1 
    clear dParams;
    dParams.plotLayout = {4,3};
    dParams.RTBins = 10;
    h = figure(1);
    set(h, 'units', 'inches', 'position', [20, 5, 8.5, 11]);
    clf;
    inited = true;
  else
    inited = false;
  end
  
  % additionally, we always force subjectNumber to be scalar
  if ~isempty(file) && ~isscalar(file.subjectNumber)
      file.subjectNumber = file.subjectNumber(end);
  end
end