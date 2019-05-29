% This function computes the average for each condition for the given subject group
function AvgData = MP_GetGroupAverage(curr_folder, subject_group, events, prefix, modeDesc)

    % Get data for selected mode
    fileList = dir(fullfile( curr_folder, strcat('*','_',modeDesc,'.mat') ));
    fileNames = sort({fileList.name});

    subj_count = size(subject_group,1);
    
    % Conditions for file selection
    c1 = contains( fileNames, 'MPData_' );
    c2 = contains( fileNames, prefix );
    c3 = contains( fileNames, modeDesc );
    
    for s=1:subj_count % Loop over subjects

        % Get data files for current subject   
        
        c4 = contains( fileNames, strcat('Subject_',subject_group(s,:)) ); 
 
        cond_num = size(events,2);
        
        for i=1:cond_num % Loop over conditions for current subject

            c5 = contains( fileNames, strcat('Condition_',num2str(events(i))) );
            fileIndices = find( c1 & c2 & c3 & c4 & c5);

            if( size(fileIndices,1) > 1 || size(fileIndices,1) == 0)
                error('Error in data file selection!');
                return;        
            end       
            
            currFile = fileNames{1,fileIndices};
            
            data = load( strcat(curr_folder,'/',currFile) );
            data = data.MPoutput;
            energy = data.meanEnergy;

            if( exist('EnergySum') == 0 ) % initialize sum matrices in first iteration
                initMat = zeros(size(energy));
                for c=1:cond_num
                    EnergySum{c} = initMat;
                end
            end

            % Accumulate sum for current condition
            EnergySum{1,i}  = EnergySum{1,i} + energy;

        end
      
    end
    
    % Calculate average for each condition
    for j = 1:size(EnergySum,2)
        EnergySum{1,j} = EnergySum{1,j}/subj_count;
    end
    
    AvgData.Avg = EnergySum;
    % Assumption: Those are the same for all subjects
    AvgData.Fs = data.Fs;
    AvgData.frequency = data.frequency;
    AvgData.time = data.time;
    
end