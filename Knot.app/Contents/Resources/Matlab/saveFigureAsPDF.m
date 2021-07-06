function saveFigureAsPDF(figureNum, fileName)

if ~ishandle(figureNum)         % first check that the figure exists
    return
end;

printParams = {
    'size', 7.*[1 0.75] ...     % [width height]
    'units', 'inches', ...      % units of size
    'resolution', 400, ...      % dpi
    'fileFormat', 'pdf', ...    % see -d options to PRINT
    'Renderer', '', ...
    'boundingBox', [], ...
    'CMYKColor', false, ...
    'PrintUI', false, ...
};

% save old paper pos / units
% note order is important here: e.g. units must come before pos

fieldsToSave = {'PaperType', 'PaperUnits', 'PaperOrientation', 'PaperPosition', 'PaperPositionMode', 'PaperSize'};
savedCell = get(figureNum, fieldsToSave);
set(figureNum, 'PaperUnits', 'inches', 'PaperOrientation', 'Portrait', 'PaperPosition', [0.5 0.5 7.5, 10], ...
    'PaperSize', [8.5, 11]);
pOptions = {'-dpdf', '-r400', '-cmyk', '-noui'};
print(figureNum, fileName, pOptions{:});
set(figureNum, fieldsToSave, savedCell);
