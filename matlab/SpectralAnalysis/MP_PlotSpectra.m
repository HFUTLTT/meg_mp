function fig = MP_PlotSpectra(SpectraList,figSpecConfig)
   
    fig = figure('Position', get(0, 'Screensize')); 
    leg = [];
    
    for i=1:size(SpectraList,1) % Loop over data list
        
        % Plot spectra
        plot_data = SpectraList(i,:);
        
        ax = gca;
        plot1 = plot (ax, plot_data.frequency, plot_data.dEnRatio); 
        xlim(figSpecConfig.freqLimsHz);
        xlabel('Frequency (Hz)');
        ylabel('Change in Power for Post-Stimulus Signal');
        hold on;
        
        set(plot1,'LineWidth',2);
        leg = [leg; strcat('Condition: ', convertCharsToStrings(plot_data.category) ) ];   
        
    end

    ax.XGrid = 'on';
    ax.XTick = [figSpecConfig.freqLimsHz(1):5:figSpecConfig.freqLimsHz(2)];
    ax.XMinorTick = 'on';        
    ax.YMinorTick = 'on';
    set(findall(gcf,'-property','FontSize'),'FontSize',figSpecConfig.fontSizeSmall);
    title(figSpecConfig.title);
    legend(leg,'Location','northeast');
    hold off;

 
end