%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR GENDER DIFFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify colours
colours.two = {[0.7 0.7 0.7]; [0.7 0.7 0.7]};
names = {'Women', 'Men'};
input_options = {'color', colours.two, 'nofig', 'names', names};

% Set the figure inputs
FDTresults.filters{:,1} = FDT.data.summary.avgFilter(FDT.GLMs.SAGI(:,2) == 1);
FDTresults.filters{:,2} = FDT.data.summary.avgFilter(FDT.GLMs.SAGI(:,2) == 0);
FDTresults.bias{:,1} = FDT.fit.SA.mratio.modelFit.c1(FDT.GLMs.SAGI(:,2) == 1);
FDTresults.bias{:,2} = FDT.fit.SA.mratio.modelFit.c1(FDT.GLMs.SAGI(:,2) == 0);
FDTresults.confidence{:,1} = FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 1);
FDTresults.confidence{:,2} = FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 0);
FDTresults.mratio{:,1} = FDT.fit.SA.mratio.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 1);
FDTresults.mratio{:,2} = FDT.fit.SA.mratio.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 0);

% Make violin figure for FDT results split by gender
figure;
    set(gcf, 'Position', [0 0 300 800])
    % Plot threshold filter sensitivity
    subplot(4,1,1);
    barplot_columns(FDTresults.filters, input_options{:});
    ylim([0, 12]);
    ylabel('Threshold filter', 'FontSize', 14);
    title('SENSITIVITY', 'FontSize', 20);
    % Plot decision bias
    subplot(4,1,2);
    barplot_columns(FDTresults.bias, input_options{:});
    ylim([-1.8, 1.8]);
    ylabel('Decision bias \it c', 'FontSize', 14);
    title('DECISION BIAS', 'FontSize', 20);
    % Plot average confidence
    subplot(4,1,3);
    barplot_columns(FDTresults.confidence, input_options{:});
    ylim([0, 12]);
    ylabel('Confidence score', 'FontSize', 14);
    title('META. BIAS', 'FontSize', 20);
    % Plot Mratio
    subplot(4,1,4);
    barplot_columns(FDTresults.mratio, input_options{:});
    ylim([0.5, 1.5]);
    ylabel('Mratio', 'FontSize', 14);
    title('META. INSIGHT', 'FontSize', 20);
    print(FDT.settings.names.figure1output, '-dtiff');
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR ANXIETY SLOPES AND GENDER INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anxiety slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.SA(:,1), FDT.data.summary.avgFilter, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
%     rf0 = refline(0, FDT.fit.SA.sensitivity.coefficients(1));
%     rf0.LineStyle = ':';
%     rf0.Color = [0.7 0.7 0.7];
%     rf0.LineWidth = 2;
    rf1 = refline(FDT.fit.SA.sensitivity.coefficients(2), FDT.fit.SA.sensitivity.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Thresh filter number'); 
    title('SENSITIVITY vs ANXIETY', 'fontsize', 20);
    % Decision bias
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.SA(:,1), FDT.fit.SA.mratio.modelFit.c1, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
%     rf0 = refline(0, FDT.fit.SA.bias.coefficients(1));
%     rf0.LineStyle = ':';
%     rf0.Color = [0.7 0.7 0.7];
%     rf0.LineWidth = 2;
    rf1 = refline(FDT.fit.SA.bias.coefficients(2), FDT.fit.SA.bias.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c'); 
    title('DECISION BIAS vs ANXIETY', 'fontsize', 20);
    % Sensibility
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.SA(:,1), FDT.data.summary.avgConfidence, 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
%     rf0 = refline(0, FDT.fit.SA.sensibility.coefficients(1));
%     rf0.LineStyle = ':';
%     rf0.Color = [0.7 0.7 0.7];
%     rf0.LineWidth = 2;
    rf1 = refline(FDT.fit.SA.sensibility.coefficients(2), FDT.fit.SA.sensibility.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['Beta = ', sprintf('%.2f', FDT.fit.SA.sensibility.coefficients(2)), '*']};
    text(3, 4, str, 'FontSize', 16, 'Color', 'r');
    title('META. BIAS vs ANXIETY', 'fontsize', 20);
    % MRatio
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 1)), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    scatter(ax4, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 0)), 25, 'filled', 'MarkerFaceColor', 'b');
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
%     rf0 = refline(0, FDT.fit.SAG.mratio.SAGIS.modelFit.mu_logMratio);
%     rf0.LineStyle = ':';
%     rf0.Color = [0.7 0.7 0.7];
%     rf0.LineWidth = 2;
    rf1 = refline((FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3), (FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta2));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    rf2 = refline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1, FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio);
    rf2.LineStyle = '--';
    rf2.Color = 'b';
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. State Anxiety score');
    str = {['Beta = ', sprintf('%.2f', (FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3)), '*']};
    text(3, -0.5, str, 'FontSize', 16, 'Color', 'r');
    title('META. INSIGHT vs ANXIETY', 'fontsize', 20);
    lgd = legend('Women','Men');
    lgd.Location = 'southwest';
    print(FDT.settings.names.figure2output, '-dtiff')

