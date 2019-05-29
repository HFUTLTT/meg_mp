% Set up MP ouput folders
function MP_SetupOuputFolder(delete, folderPath)

    if (delete == 1 && exist(folderPath,'file') > 0) % Clear if it exists and delete==1
      rmdir(folderPath, 's');    
    end
    
    if ( exist(folderPath,'file') == 0 ) % Create output folder if it does not already exist
        mkdir(folderPath);
    end    
    
end