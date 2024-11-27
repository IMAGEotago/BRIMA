%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR STATE ANXIETY SLOPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove missing values
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));
Dep = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
end
[idxNaN_row_D col] = find(isnan(Dep));
sensitivity = FDT.data.summary.avgFilter;
sensitivity(:,idxNaN_row_D) = [];
confidence = FDT.data.summary.avgConfidence;
confidence(:,idxNaN_row_D) = [];


% Depression slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.D(:,1), sensitivity, 25, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.D.sensitivity.coefficients(2), FDT.fit.D.sensitivity.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Threshold filter'); 
    title('SENSITIVITY', 'fontsize', 18);
    % Decision bias
    ax2 = subplot(4,1,2);
    scatter(ax2, FDT.GLMs.D(:,1), FDT.fit.D.mratio.modelFit.c1, 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.D.bias.coefficients(2), FDT.fit.D.bias.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'Decision bias  \it c');
    str = {['Beta = ', sprintf('%.2f', FDT.fit.D.bias.coefficients(2)), '*']};
    text(3, 0.8, str, 'FontSize', 12, 'Color', 'r');
    title('DECISION BIAS', 'fontsize', 18);
    % Metacognitive bias
    ax3 = subplot(4,1,3);
    scatter(ax3, FDT.GLMs.D(:,1), confidence, 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
    rf1 = refline(FDT.fit.D.sensibility.coefficients(2), FDT.fit.D.sensibility.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['Beta = ', sprintf('%.2f', FDT.fit.D.sensibility.coefficients(2)), '*']};
    text(3, 9, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    % MRatio
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.D(:,1), log(FDT.fit.D.mratio.modelFit.Mratio), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
    rf1 = refline(FDT.fit.D.mratio.modelFit.mu_beta1, FDT.fit.D.mratio.modelFit.mu_logMratio);
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. Deoression score');
    str = {['Beta = ', sprintf('%.2f', FDT.fit.D.mratio.modelFit.mu_beta1), '#']};
    text(3, -0.75, str, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    print(FDT.settings.names.figure.D1, '-dtiff');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR DEPRESSION GENDER DIFFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify colours
colours.two = {[0.7 0.7 0.7]; [0.7 0.7 0.7]};
names = {'Women', 'Men'};
input_options = {'color', colours.two, 'nofig', 'names', names};

% Set the figure inputs
FDTresults.filters{:,1} = FDT.data.summary.avgFilter(FDT.GLMs.DGI(:,2) == 1);
FDTresults.filters{:,2} = FDT.data.summary.avgFilter(FDT.GLMs.DGI(:,2) == 0);
FDTresults.bias{:,1} = FDT.fit.D.mratio.modelFit.c1(FDT.GLMs.DGI(:,2) == 1);
FDTresults.bias{:,2} = FDT.fit.D.mratio.modelFit.c1(FDT.GLMs.DGI(:,2) == 0);
FDTresults.confidence{:,1} = FDT.data.summary.avgConfidence(FDT.GLMs.DGI(:,2) == 1);
FDTresults.confidence{:,2} = FDT.data.summary.avgConfidence(FDT.GLMs.DGI(:,2) == 0);
FDTresults.mratio{:,1} = FDT.fit.D.mratio.modelFit.Mratio(FDT.GLMs.DGI(:,2) == 1);
FDTresults.mratio{:,2} = FDT.fit.D.mratio.modelFit.Mratio(FDT.GLMs.DGI(:,2) == 0);

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
    print(FDT.settings.names.figure.D2, '-dtiff');
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR DEPRESSION SLOPES AND GENDER INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Depression slopes
figure
    set(gcf, 'Position', [0 0 500 800]);
    % Sensitivity: gender
    ax1 = subplot(4,1,1);
    scatter(ax1, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 0),1), FDT.data.summary.avgFilter(FDT.GLMs.DGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax1, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 1),1), FDT.data.summary.avgFilter(FDT.GLMs.DGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [0 10])
    rf1 = refline(FDT.fit.DG.sensitivity.DGI.coefficients(2), FDT.fit.DG.sensitivity.DGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.DG.sensitivity.DGI.coefficients(2) + FDT.fit.DG.sensitivity.DGI.coefficients(4)), (FDT.fit.DG.sensitivity.DGI.coefficients(1) + FDT.fit.DG.sensitivity.DGI.coefficients(3)));
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
    scatter(ax2, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 0),1), FDT.fit.DG.mratio.DGI.modelFit.c1(FDT.GLMs.DGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax2, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 1),1), FDT.fit.DG.mratio.DGI.modelFit.c1(FDT.GLMs.DGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1.5 1.5])
    rf1 = refline(FDT.fit.DG.bias.DGI.coefficients(2), FDT.fit.DG.bias.DGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = [0.7 0.7 0.7];
    rf1.LineWidth = 2;
    rf2 = refline((FDT.fit.DG.bias.DGI.coefficients(2) + FDT.fit.DG.bias.DGI.coefficients(4)), (FDT.fit.DG.bias.DGI.coefficients(1) + FDT.fit.DG.bias.DGI.coefficients(3)));
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
    scatter(ax3, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 1),1), FDT.data.summary.avgConfidence(FDT.GLMs.DGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax3, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 0),1), FDT.data.summary.avgConfidence(FDT.GLMs.DGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax3, [-2 5])
    ylim(ax3, [2 11])
    rf2 = refline((FDT.fit.DG.sensibility.DGI.coefficients(2) + FDT.fit.DG.sensibility.DGI.coefficients(4)), (FDT.fit.DG.sensibility.DGI.coefficients(1) + FDT.fit.DG.sensibility.DGI.coefficients(3)));
    rf2.LineStyle = '--';
    rf2.Color = [0.7 0.7 0.7];
    rf2.LineWidth = 2;
    rf1 = refline(FDT.fit.DG.sensibility.DGI.coefficients(2), FDT.fit.DG.sensibility.DGI.coefficients(1));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax3,'YLabel'), 'String', 'Avg. Confidence'); 
    str = {['M Beta = ', sprintf('%.2f', (FDT.fit.DG.sensibility.DGI.coefficients(2))), '*']};
    text(3, 3.2, str, 'FontSize', 12, 'Color', 'r');
    title('METACOGNITIVE BIAS', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'northeast';
    % MRatio: gender
    ax4 = subplot(4,1,4);
    scatter(ax4, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 0),1), log(FDT.fit.DG.mratio.DGI.modelFit.Mratio(FDT.GLMs.DGI(:,2) == 0)), 25, 'filled', 'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on
    scatter(ax4, FDT.GLMs.DGI((FDT.GLMs.DGI(:,2) == 1),1), log(FDT.fit.DG.mratio.DGI.modelFit.Mratio(FDT.GLMs.DGI(:,2) == 1)), 25, 'filled', 'MarkerFaceColor', 'r');
    box on
    xlim(ax4, [-2 5])
    ylim(ax4, [-1 0.3])
    rf1 = refline((FDT.fit.DG.mratio.DGI.modelFit.mu_beta1 + FDT.fit.DG.mratio.DGI.modelFit.mu_beta3), (FDT.fit.DG.mratio.DGI.modelFit.mu_logMratio + FDT.fit.DG.mratio.DGI.modelFit.mu_beta2));
    rf1.LineStyle = '--';
    rf1.Color = 'r';
    rf1.LineWidth = 2;
    rf2 = refline(FDT.fit.DG.mratio.DGI.modelFit.mu_beta1, FDT.fit.DG.mratio.DGI.modelFit.mu_logMratio);
    rf2.LineStyle = '--';
    rf2.Color = [0.7 0.7 0.7];
    rf2.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax4,'YLabel'), 'String', 'logMratio'); 
    set(get(ax4,'XLabel'), 'String', 'Std. Depression score');
    str1 = {['W Beta = ', sprintf('%.2f', (FDT.fit.DG.mratio.DGI.modelFit.mu_beta1 + FDT.fit.DG.mratio.DGI.modelFit.mu_beta3)), '#']};
    text(-1.7, -0.9, str1, 'FontSize', 12, 'Color', 'r');
    title('INSIGHT', 'fontsize', 18);
    lgd = legend('Men','Women');
    lgd.Location = 'northeast';
    print(FDT.settings.names.figure.D3, '-dtiff');

