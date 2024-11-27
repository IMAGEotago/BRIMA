%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FIGURE FOR STATE ANXIETY SLOPES AND GENDER INTERACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Anxiety slopes for women
figure
    set(gcf, 'Position', [0 0 350 400]);
    % Metacognitive bias: gender
    ax1 = subplot(2,1,1);
    scatter(ax1, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 1), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [2 11])
    rf = refline((FDT.fit.SAG.sensibility.SAGI.coefficients(2) + FDT.fit.SAG.sensibility.SAGI.coefficients(4)), (FDT.fit.SAG.sensibility.SAGI.coefficients(1) + FDT.fit.SAG.sensibility.SAGI.coefficients(3)));
    rf.LineStyle = '--';
    rf.Color = 'r';
    rf.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Avg. Confidence');
    set(get(ax1,'XLabel'), 'String', 'Std. State Anxiety score');
    
    % MRatio: gender
    ax2 = subplot(2,1,2);
    scatter(ax2, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 1),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 1)), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1 0.3])
    rf = refline((FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3), (FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta2));
    rf.LineStyle = '--';
    rf.Color = 'r';
    rf.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'logMratio'); 
    set(get(ax2,'XLabel'), 'String', 'Std. State Anxiety score');
    print('graphicalAbstract-figure-women', '-djpeg')

% Anxiety slopes for men
figure
    set(gcf, 'Position', [0 0 350 400]);
    % Metacognitive bias: gender
    ax1 = subplot(2,1,1);
    scatter(ax1, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), FDT.data.summary.avgConfidence(FDT.GLMs.SAGI(:,2) == 0), 25, 'filled', 'MarkerFaceColor', 'r');
    hold on
    box on
    xlim(ax1, [-2 5])
    ylim(ax1, [2 11])
    rf = refline(FDT.fit.SAG.sensibility.SAGI.coefficients(2), FDT.fit.SAG.sensibility.SAGI.coefficients(1));
    rf.LineStyle = '--';
    rf.Color = 'r';
    rf.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax1,'YLabel'), 'String', 'Avg. Confidence'); 
    set(get(ax1,'XLabel'), 'String', 'Std. State Anxiety score');

    % MRatio: gender
    ax2 = subplot(2,1,2);
    scatter(ax2, FDT.GLMs.SAGI((FDT.GLMs.SAGI(:,2) == 0),1), log(FDT.fit.SAG.mratio.SAGI.modelFit.Mratio(FDT.GLMs.SAGI(:,2) == 0)), 25, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
    hold on
    box on
    xlim(ax2, [-2 5])
    ylim(ax2, [-1 0.3])
    rf = refline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1, FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio);
    rf.LineStyle = '--';
    rf.Color = [0.3 0.3 0.3];
    rf.LineWidth = 2;
    set(gca, 'FontSize', 14)
    set(get(ax2,'YLabel'), 'String', 'logMratio'); 
    set(get(ax2,'XLabel'), 'String', 'Std. State Anxiety score');
    print('graphicalAbstract-figure-men', '-djpeg')