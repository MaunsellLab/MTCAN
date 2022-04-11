

fileName = 'testing_220325001.ns5';                         % .ns5 file 
cd '/Users/chery/Documents/Grad School/Maunsell Lab/Analysis/MTCAN/NEV/';  % directory for .nev file
openNSx('read', fileName, 'uv');
myData = NS5.Data;

saveName = strrep(fileName, '.ns5', '.bin');
fileID = fopen(saveName, 'w');
fwrite(fileID, myData, 'int16')
fclose(fileID);
