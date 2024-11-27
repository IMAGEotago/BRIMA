%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR STATE ANXIETY SLOPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% State anxiety slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.SA(:,1), FDT.data.summary.avgFilter, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.SA.sensitivity.coefficients(2), FDT.fit.SA.sensitivity.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Threshold filter'); 
    title('SENSITIVITY', 'fontsize', 18);
    % Decision bias
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.SA(:,1), FDT.fit.SA.mratio.modelFit.c1, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.SA.bias.coefficients(2), FDT.fit.SA.bias.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c'); 
    title('DECISION BIAS', 'fontsize', 18);
    % Metacognitive bias
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.SA(:,1), FDT.data.summary.avgConfidence, 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
    rf1 = refline(FDT.fit.SA.sensibility.coefficients(2), FDT.fit.SA.sensibility.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['Beta = ', sprintf('%.2f', FDT.fit.SA.sensibility.coefficients(2)), '*']};
    text(3, 4, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    % MRatio
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.SA(:,1), log(FDT.fit.SA.mratio.modelFit.Mratio), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
    rf1 = refline(FDT.fit.SA.mratio.modelFit.mu_beta1, FDT.fit.SA.mratio.modelFit.mu_logMratio);
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. State Anxiety score');
    str = {['Beta = ', sprintf('%.2f', FDT.fit.SA.mratio.modelFit.mu_beta1), '#']};
    text(3, -0.75, str, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    print(FDT.settings.names.figure.SA1, '-dtiff')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR STATE ANXIETY GENDER DIFFERENCES
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
    title('SENSITIVITY', 'FontSize', 18);
    % Plot decision bias
    subplot(4,1,2);
    barplot_columns(FDTresults.bias, input_options{:});
    ylim([-1.8, 1.8]);
    ylabel('Decision bias \it c', 'FontSize', 14);
    title('DECISION BIAS', 'FontSize', 18);
    % Plot average confidence
    subplot(4,1,3);
    barplot_columns(FDTresults.confidence, input_options{:});
    ylim([0, 12]);
    ylabel('Confidence score', 'FontSize', 14);
    title('META. BIAS', 'FontSize', 18);
    % Plot Mratio
    subplot(4,1,4);
    barplot_columns(FDTresults.mratio, input_options{:});
    ylim([0.5, 1.5]);
    ylabel('Mratio', 'FontSize', 14);
    title('INSIGHT', 'FontSize', 18);
    print(FDT.settings.names.figure.SA2, '-dtiff');
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR STATE ANXIETY SLOPES AND GENDER INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anxiety slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity: gender
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), FDT.data.summary.avgFilter(FDT.GLMs.SAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax1, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), FDT.data.summary.avgFilter(FDT.GLMs.SAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.SAG.sensitivity.SAGI.coefficients(2), FDT.fit.SAG.sensitivity.SAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.SAG.sensitivity.SAGI.coefficients(2) + FDT.fit.SAG.sensitivity.SAGI.coefficients(4)), (FDT.fit.SAG.sensitivity.SAGI.coefficients(1) + FDT.fit.SAG.sensitivity.SAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.3 0.3 0.3];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Threshold filter'); 
    title('SENSITIVITY', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'east';
    % Decision bias: gender
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), FDT.fit.SAG.mratio.SAGI.modelFit.c1(FDT.GLMs.SAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax2, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), FDT.fit.SAG.mratio.SAGI.modelFit.c1(FDT.GLMs.SAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.SAG.bias.SAGI.coefficients(2), FDT.fit.SAG.bias.SAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.SAG.bias.SAGI.coefficients(2) + FDT.fit.SAG.bias.SAGI.coefficients(4)), (FDT.fit.SAG.bias.SAGI.coefficients(1) + FDT.fit.SAG.bias.SAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.3 0.3 0.3];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c'); 
    title('DECISION BIAS', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'southeast';
    % Metacognitive bias: gender
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax3, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [-2 11])
    rf1 = refline(FDT.fit.SAG.sensibility.SAGI.coefficients(2), FDT.fit.SAG.sensibility.SAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.SAG.sensibility.SAGI.coefficients(2) + FDT.fit.SAG.sensibility.SAGI.coefficients(4)), (FDT.fit.SAG.sensibility.SAGI.coefficients(1) + FDT.fit.SAG.sensibility.SAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = 'r';
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['W Beta = ', sprintf('%.2f', (FDT.fit.SAG.sensibility.SAGI.coefficients(2) + FDT.fit.SAG.sensibility.SAGI.coefficients(4))), '*']};
    text(3, 2, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'southwest';
    % MRatio: gender
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 0)), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax4, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 1)), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
    rf1 = refline((FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3), (FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta2));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    rf2 = refline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1, FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio);
    rf2.LineStyle = '--';
    rf2.Color = [0.7 0.7 0.7];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. State Anxiety score');
    str1 = {['W Beta = ', sprintf('%.2f', (FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3)), '#']};
    text(3, -0.3, str1, 'FontSize', 12, 'Color', 'r');
    str2 = {['W>M Beta = ', sprintf('%.2f', (FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3)), '#']};
    text(3, -0.5, str2, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'southwest';
    print(FDT.settings.names.figure.SA3, '-dtiff')

