%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR TRAIT ANXIETY SLOPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove missing values
idxTA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-T'));
TA = [];
for study = 1:4
    TA = [TA; FDT.data.raw.extra{study}(:,idxTA)];
end
[idxNaN_row_TA col] = find(isnan(TA));
sensitivity = FDT.data.summary.avgFilter;
sensitivity(:,idxNaN_row_TA) = [];
confidence = FDT.data.summary.avgConfidence;
confidence(:,idxNaN_row_TA) = [];

% Trait anxiety slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.TA(:,1), sensitivity, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.TA.sensitivity.coefficients(2), FDT.fit.TA.sensitivity.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Threshold filter'); 
    title('SENSITIVITY', 'fontsize', 18);
    % Decision bias
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.TA(:,1), FDT.fit.TA.mratio.modelFit.c1, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.TA.bias.coefficients(2), FDT.fit.TA.bias.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c'); 
    title('DECISION BIAS', 'fontsize', 18);
    % Metacognitive bias
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.TA(:,1), confidence, 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
    rf1 = refline(FDT.fit.TA.sensibility.coefficients(2), FDT.fit.TA.sensibility.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['Beta = ', sprintf('%.2f', FDT.fit.TA.sensibility.coefficients(2)), '*']};
    text(3, 3.4, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    % MRatio
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.TA(:,1), log(FDT.fit.TA.mratio.modelFit.Mratio), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
    rf1 = refline(FDT.fit.TA.mratio.modelFit.mu_beta1, FDT.fit.TA.mratio.modelFit.mu_logMratio);
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. Trait Anxiety score');
    str = {['Beta = ', sprintf('%.2f', FDT.fit.TA.mratio.modelFit.mu_beta1), '#']};
    text(3, -0.75, str, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    print(FDT.settings.names.figure.TA1, '-dtiff')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR TRAIT ANXIETY GENDER DIFFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify colours
colours.two = {[0.7 0.7 0.7]; [0.7 0.7 0.7]};
names = {'Women', 'Men'};
input_options = {'color', colours.two, 'nofig', 'names', names};

% Set the figure inputs
FDTresults.filters{:,1} = FDT.data.summary.avgFilter(FDT.GLMs.TAGI(:,2) == 1);
FDTresults.filters{:,2} = FDT.data.summary.avgFilter(FDT.GLMs.TAGI(:,2) == 0);
FDTresults.bias{:,1} = FDT.fit.TA.mratio.modelFit.c1(FDT.GLMs.TAGI(:,2) == 1);
FDTresults.bias{:,2} = FDT.fit.TA.mratio.modelFit.c1(FDT.GLMs.TAGI(:,2) == 0);
FDTresults.confidence{:,1} = FDT.data.summary.avgConfidence(FDT.GLMs.TAGI(:,2) == 1);
FDTresults.confidence{:,2} = FDT.data.summary.avgConfidence(FDT.GLMs.TAGI(:,2) == 0);
FDTresults.mratio{:,1} = FDT.fit.TA.mratio.modelFit.Mratio(FDT.GLMs.TAGI(:,2) == 1);
FDTresults.mratio{:,2} = FDT.fit.TA.mratio.modelFit.Mratio(FDT.GLMs.TAGI(:,2) == 0);

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
    print(FDT.settings.names.figure.TA2, '-dtiff');
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR TRAIT ANXIETY SLOPES AND GENDER INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trait anxiety slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity: gender
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 0),1), FDT.data.summary.avgFilter(FDT.GLMs.TAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax1, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 1),1), FDT.data.summary.avgFilter(FDT.GLMs.TAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.TAG.sensitivity.TAGI.coefficients(2), FDT.fit.TAG.sensitivity.TAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.TAG.sensitivity.TAGI.coefficients(2) + FDT.fit.TAG.sensitivity.TAGI.coefficients(4)), (FDT.fit.TAG.sensitivity.TAGI.coefficients(1) + FDT.fit.TAG.sensitivity.TAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.3 0.3 0.3];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Threshold filter'); 
    title('SENSITIVITY', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'northeast';
    % Decision bias: gender
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 0),1), FDT.fit.TAG.mratio.TAGI.modelFit.c1(FDT.GLMs.TAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax2, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 1),1), FDT.fit.TAG.mratio.TAGI.modelFit.c1(FDT.GLMs.TAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.TAG.bias.TAGI.coefficients(2), FDT.fit.TAG.bias.TAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.TAG.bias.TAGI.coefficients(2) + FDT.fit.TAG.bias.TAGI.coefficients(4)), (FDT.fit.TAG.bias.TAGI.coefficients(1) + FDT.fit.TAG.bias.TAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.3 0.3 0.3];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c'); 
    title('DECISION BIAS', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'northeast';
    % Metacognitive bias: gender
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 1),1), FDT.data.summary.avgConfidence(FDT.GLMs.TAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax3, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 0),1), FDT.data.summary.avgConfidence(FDT.GLMs.TAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
    rf2 = refline((FDT.fit.TAG.sensibility.TAGI.coefficients(2) + FDT.fit.TAG.sensibility.TAGI.coefficients(4)), (FDT.fit.TAG.sensibility.TAGI.coefficients(1) + FDT.fit.TAG.sensibility.TAGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.7 0.7 0.7];
    rf2.LineWidth = 2;
    rf1 = refline(FDT.fit.TAG.sensibility.TAGI.coefficients(2), FDT.fit.TAG.sensibility.TAGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['M Beta = ', sprintf('%.2f', (FDT.fit.TAG.sensibility.TAGI.coefficients(2))), '*']};
    text(3, 3.5, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    lgd = legend('Women','Men');
    lgd.Location = 'northeast';
    % MRatio: gender
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 0),1), log(FDT.fit.TAG.mratio.TAGI.modelFit.Mratio(FDT.GLMs.TAGI(:,2) == 0)), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax4, FDT.GLMs.TAGI((FDT.GLMs.TAGI(:,2) == 1),1), log(FDT.fit.TAG.mratio.TAGI.modelFit.Mratio(FDT.GLMs.TAGI(:,2) == 1)), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-0.6 0.3])
    rf1 = refline((FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1 + FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3), (FDT.fit.TAG.mratio.TAGI.modelFit.mu_logMratio + FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta2));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    rf2 = refline(FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1, FDT.fit.TAG.mratio.TAGI.modelFit.mu_logMratio);
    rf2.LineStyle = '--';
    rf2.Color = [0.7 0.7 0.7];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. Trait Anxiety score');
    str = {['W Beta = ', sprintf('%.2f', (FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1 + FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3)), '#']};
    text(-1.7, -0.5, str, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'northeast';
    print(FDT.settings.names.figure.TA3, '-dtiff')

