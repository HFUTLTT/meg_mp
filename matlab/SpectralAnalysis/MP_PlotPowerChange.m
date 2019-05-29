function figPC = MP_PlotPowerChange(dEnDB, figConfig)
   
    tfData{1} = dEnDB; % Change in power

    hTF=cell(1,1);

    figPC = figure('Position', get(0, 'Screensize'));
    for i=1:1
        hTF{i} = subplot(1,1,1);  % Figure handle
        imagesc(figConfig.time, figConfig.frequency, tfData{i}, 'Parent', hTF{i});
        set(hTF{i},'YDir','normal');
        axis(hTF{i},[figConfig.timeLimsS figConfig.freqLimsHz]);
        if i==1
            caxis(hTF{i},figConfig.cLimsDiff);
        else
            caxis(hTF{i},figConfig.cLims);
        end
    end

    % Figure Details
    titleList =[{'Change in Power', '', ''}];
    plotNumberList = {'A','B','C'};
    for i=1:1
        title(hTF{i},titleList{i},'FontSize',figConfig.fontSizeSmall); % Titles
        % Labels and Fonts
        set(hTF{i},'XTick',[figConfig.time(1):0.1:figConfig.time(end)],'YTick',(0:0.1:1)*figConfig.freqLimsHz(2),'FontSize',figConfig.fontSizeSmall);
        xlabel(hTF{i}, 'Time (s)', 'FontSize', figConfig.fontSizeSmall);
        box(hTF{i}, 'on'); % Boxes
        colorbar('peer', hTF{i}, 'location', 'northoutside'); % Colorbars
        text(-0.1,1.5, plotNumberList{i}, 'unit','normalized','fontsize',figConfig.fontSizeLarge,'Parent',hTF{i}); % Text
    end
    ylabel(hTF{1},'Frequency (Hz)','FontSize',figConfig.fontSizeSmall);
    
    title(figConfig.title);
    
end