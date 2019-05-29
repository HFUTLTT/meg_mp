%% Environment Setup

% Close and clear everything
clear all;
close all;
clc;

%% Path and Folder Setup
curr_folder = "convert/";
SourceFileType = ".fig";
TargetFileType = "mpeg";
% Get data for selected mode
fileList = dir(fullfile( curr_folder, strcat('*',SourceFileType) ));
fileNames = sort({fileList.name});

%%
n = size(fileNames,2);

for i=1:n
    
    fname = fileNames(i);
    file2load = strcat(curr_folder,fname);
    %data = load( file2load );  % for .mat files
    fig = openfig( file2load, 'new', 'invisible'); % for .fig file
    outputName = erase(fname,SourceFileType);
    
    % for .mat files
    %{
    if(isfield(data,'fig1'))
        saveas(data.fig1, strcat(curr_folder, outputName,'_fig1'), 'jpeg');
    end
    
    if(isfield(data,'fig2'))
        saveas(data.fig2, strcat(curr_folder, outputName,'_fig2'), 'jpeg');
    end
    
    if(isfield(data,'fig3'))
        disp('Something is wrong!')
    end
    %}
    
    saveas(fig, strcat(curr_folder, outputName), 'jpeg');
    
end

%% 
disp('Done!')
%close all;