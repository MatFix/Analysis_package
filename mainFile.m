function renameFiles(folder,PythonScript)
%% MATLAB wrapper for a Python function
%
%  Calls the Python script batchRenamer.py (Requires Python 3.5 or above)
%

folderRename = ['"' folder];

commandStr = ['python ' '"' pwd '\' PythonScript '" -r ' folderRename];
[status] = system(commandStr);

if status == 0
    disp('Renaming successful')
else 
    error('Error: could not rename')
end

end