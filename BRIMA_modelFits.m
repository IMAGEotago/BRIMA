%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BRIMA DATA FDT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Olivia Harrison
% Created: 12/04/2022
% -------------------------------------------------------------------------
% TO RUN:                 BRIMA_modelFits
% OUTPUTS:              = The analysis of the empirical FDT datasets
% -------------------------------------------------------------------------
% DESCRIPTION:
% This analysis performs the empirical data assessment of the FDT from a
% combined dataset of individuals collected in Bath, Birmingham, Zurich and
% Oxford.
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FDT = BRIMA_modelFits

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save specified options
FDT.settings = BRIMA_setOptions;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT THE RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the datasheets for each study
for study = 1:4
    [FDT.data.raw.extra{study}, tempIDs{study}] = xlsread(FDT.settings.names.dataFile, FDT.settings.names.sheetNames{study});
    FDT.data.IDs{study} = tempIDs{study}(2:end,1);
end

% Load the FDT data from Bath, Birmingham and Zurich
for study = 1:3
    for input = 1:length(FDT.data.IDs{study})
        if strcmp(FDT.settings.names.sheetNames{study}, 'BATH')
            ID = FDT.data.IDs{study}{input}(end-2:end);
        elseif strcmp(FDT.settings.names.sheetNames{study}, 'BIRMINGHAM')
            ID = [FDT.data.IDs{study}{input}, 'V1'];
        elseif strcmp(FDT.settings.names.sheetNames{study}, 'ZURICH')
            ID = FDT.data.IDs{study}{input}(end-3:end);
        end
        FDT.data.raw.FDT{study}{input} = load(fullfile(FDT.settings.paths.FDTdata, FDT.settings.names.sheetNames{study}, ['filter_task_results_', ID, '.mat']));
    end
end
for input = 1:length(FDT.data.IDs{4})
    ID = FDT.data.IDs{4}{input};
    FDT.data.raw.FDT{4}{input} = xlsread(fullfile(FDT.settings.paths.FDTdata, FDT.settings.names.sheetNames{4}, [ID, '_Interoception.xlsx']));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Format the FDT data from Bath, Birmingham and Zurich
for study = 1:3
    for input = 1:length(FDT.data.IDs{study})
        % FDT.data.processed{study}{input}.filterNum = FDT.data.raw.FDT{study}{input}.results.allTrials.intensity;
        FDT.data.processed{study}{input}.response = FDT.data.raw.FDT{study}{input}.results.thresholdTrials.response;
        FDT.data.processed{study}{input}.filterNum(1:length(FDT.data.processed{study}{input}.response)) = FDT.data.raw.FDT{study}{input}.results.thresholdTrials.filterNum;
        FDT.data.processed{study}{input}.stimID = FDT.data.raw.FDT{study}{input}.results.thresholdTrials.filters(1:length(FDT.data.processed{study}{input}.filterNum));
        FDT.data.processed{study}{input}.rating = FDT.data.raw.FDT{study}{input}.results.thresholdTrials.confidence;
    end
end

% Format and bin the Oxford data
for input = 1:length(FDT.data.IDs{4})
    % Find the threshold filter number
    tempFilterNum = FDT.data.raw.FDT{4}{input}(end,1);
    tempData = FDT.data.raw.FDT{4}{input}((FDT.data.raw.FDT{4}{input}(:,1) == tempFilterNum),:);
    FDT.data.processed{4}{input}.filterNum = tempData(:,1)';
    FDT.data.processed{4}{input}.stimID = double(tempData(:,2) > 0)';
    FDT.data.processed{4}{input}.response = double(tempData(:,3) > 0)';
    FDT.data.processed{4}{input}.rating = discretize(tempData(:,4), FDT.settings.confidenceBinEdges, 'IncludedEdge', 'right')';
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE ACCURACY, NUMBER OF TRIALS, AVERAGE FILTER AND CONFIDENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the accuracy and number of trials for each participant
count = 1;
for study = 1:4
    for input = 1:length(FDT.data.IDs{study})
        FDT.data.summary.accuracy(count) = sum(FDT.data.processed{study}{input}.stimID == FDT.data.processed{study}{input}.response) / length(FDT.data.processed{study}{input}.stimID);
        FDT.data.summary.trials(count) = length(FDT.data.processed{study}{input}.stimID);
        FDT.data.summary.avgFilter(count) = mean(FDT.data.processed{study}{input}.filterNum);
        FDT.data.summary.avgConfidence(count) = mean(FDT.data.processed{study}{input}.rating);
        count = count + 1;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN TRIALS TO COUNTS SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert all data into counts
count = 1;
for study = 1:4
    for input = 1:length(FDT.data.IDs{study})
        [FDT.data.trials2counts.nR_S1{count}, FDT.data.trials2counts.nR_S2{count}] = trials2counts(FDT.data.processed{study}{input}.stimID, FDT.data.processed{study}{input}.response, FDT.data.processed{study}{input}.rating, FDT.settings.confidenceBins, 0);
        count = count + 1;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE SA GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the state anxiety variable
idxSA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-S'));

% Create the GLM for state anxiety (with regressors for study location)
SA = [];
for study = 1:4
    SA = [SA; FDT.data.raw.extra{study}(:,idxSA)];
end
FDT.GLMs.SA = [SA, zeros(length(SA),3)];
FDT.GLMs.SA(1:length(FDT.data.IDs{1}),2:4) = -1;
FDT.GLMs.SA(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),2) = 1;
FDT.GLMs.SA(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),3) = 1;
FDT.GLMs.SA(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),4) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_SA col] = find(isnan(SA));
FDT.GLMs.SA(idxNaN_row_SA,:) = [];
FDT.data.trials2counts.SA.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.SA.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.SA.nR_S1(idxNaN_row_SA) = [];
FDT.data.trials2counts.SA.nR_S2(idxNaN_row_SA) = [];

% Zscore variables
FDT.GLMs.SA(:,1) = zscore(FDT.GLMs.SA(:,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE SA MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio
FDT.fit.SA.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SA.nR_S1, FDT.data.trials2counts.SA.nR_S2, FDT.GLMs.SA', FDT.settings.fit.mcmc_params);
FDT.fit.SA.mratio.HDI.anx_corr = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.SA.mratio.HDI.anx = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SA.mratio.HDI.intercept = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
% Plot the HDI for checking
figure;
set(gcf, 'Position', [200 200 400 300])
histogram(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.SA.mratio.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SA.mratio.HDI.anx(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SA.mratio.HDI.anx(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('State anxiety Beta HDI');
set(gca, 'FontSize', 14)

% Fit the linear model for sensitivity (filter number)
filters_SA = FDT.data.summary.avgFilter;
filters_SA(idxNaN_row_SA) = [];
FDT.fit.SA.sensitivity.modelFit = fitlm(FDT.GLMs.SA, filters_SA);
FDT.fit.SA.sensitivity.coefficients = FDT.fit.SA.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.SA.sensitivity.pvalues = FDT.fit.SA.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_SA = FDT.data.summary.avgConfidence;
confidence_SA(idxNaN_row_SA) = [];
FDT.fit.SA.sensibility.modelFit = fitlm(FDT.GLMs.SA, confidence_SA);
FDT.fit.SA.sensibility.coefficients = FDT.fit.SA.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.SA.sensibility.pvalues = FDT.fit.SA.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.SA.bias.modelFit = fitlm(FDT.GLMs.SA, FDT.fit.SA.mratio.modelFit.c1);
FDT.fit.SA.bias.coefficients = FDT.fit.SA.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.SA.bias.pvalues = FDT.fit.SA.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.SA.dPrime.modelFit = fitlm(FDT.GLMs.SA, FDT.fit.SA.mratio.modelFit.d1);
FDT.fit.SA.dPrime.coefficients = FDT.fit.SA.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.SA.dPrime.pvalues = FDT.fit.SA.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE SAG GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the state anxiety variable
idxSA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-S'));

% Find the gender variable
idxGen = find(strcmp(FDT.settings.variables.summarySheet, 'GENDER'));

% Create the gender variable and GLM with interaction specified in terms matrix
SA = [];
Gender = [];
for study = 1:4
    SA = [SA; FDT.data.raw.extra{study}(:,idxSA)];
    Gender = [Gender; FDT.data.raw.extra{study}(:,idxGen)];
end
FDT.GLMs.SAG(:,1) = SA;
FDT.GLMs.SAG(:,2) = Gender;
FDT.GLMs.SAG(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_SAG col] = find(isnan(FDT.GLMs.SAG));
FDT.GLMs.SAG(idxNaN_row_SAG,:) = [];
FDT.data.trials2counts.SAG.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.SAG.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.SAG.nR_S1(idxNaN_row_SAG) = [];
FDT.data.trials2counts.SAG.nR_S2(idxNaN_row_SAG) = [];

% Zscore variables
FDT.GLMs.SAG(:,1) = zscore(FDT.GLMs.SAG(:,1));

% Create the SAG matrix with an explicit interaction term
FDT.GLMs.SAGI = FDT.GLMs.SAG;
FDT.GLMs.SAGI(:,4:6) = FDT.GLMs.SAG(:,3:5);
FDT.GLMs.SAGI(:,3) = FDT.GLMs.SAG(:,1) .* FDT.GLMs.SAG(:,2);

% Create the SAG matrix with a split interaction term
FDT.GLMs.SAGIS(:,1) = FDT.GLMs.SAG(:,1);
FDT.GLMs.SAGIS(:,2) = FDT.GLMs.SAG(:,1);
FDT.GLMs.SAGIS(:,3) = FDT.GLMs.SAG(:,2);
FDT.GLMs.SAGIS((FDT.GLMs.SAG(:,2) == 0),1) = 0;
FDT.GLMs.SAGIS((FDT.GLMs.SAG(:,2) == 1),2) = 0;
FDT.GLMs.SAGIS(:,4:6) = FDT.GLMs.SAG(:,3:5);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE SAG MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio (specified interaction term)
FDT.fit.SAG.mratio.SAGI.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAG.nR_S1, FDT.data.trials2counts.SAG.nR_S2, FDT.GLMs.SAGI', FDT.settings.fit.mcmc_params);
FDT.fit.SAG.mratio.SAGI.HDI.anxF = calc_HDI((FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:)), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.anxM = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.gender = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.interaction = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.intercept = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.SAG.mratio.SAGI.estimates.anxF = FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3;
FDT.fit.SAG.mratio.SAGI.estimates.anxM = FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1;
FDT.fit.SAG.mratio.SAGI.estimates.gender = FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta2;
FDT.fit.SAG.mratio.SAGI.estimates.interaction = FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3;
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 800 300])
subplot(1,2,1)
histogram(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_logMratio(:));
hold on
xline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_logMratio, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.intercept(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.intercept(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Intercept (Male) HDI');
set(gca, 'FontSize', 14)
subplot(1,2,2)
histogram(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:));
hold on
xline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta2, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.gender(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.gender(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Gender (F>M) Beta HDI');
set(gca, 'FontSize', 14)
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 1200 300])
subplot(1,3,1)
histogram((FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:)));
hold on
xline((FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3), 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.anxF(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.anxF(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Female State anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,2)
histogram(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.anxM(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.anxM(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Male State anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,3)
histogram(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:));
hold on
xline(FDT.fit.SAG.mratio.SAGI.modelFit.mu_beta3, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.interaction(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG.mratio.SAGI.HDI.interaction(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Interaction Beta HDI');
set(gca, 'FontSize', 14)

% Fit the HMeta-d model for mratio (split interaction term)
FDT.fit.SAG.mratio.SAGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAG.nR_S1, FDT.data.trials2counts.SAG.nR_S2, FDT.GLMs.SAGIS', FDT.settings.fit.mcmc_params);
FDT.fit.SAG.mratio.SAGIS.HDI.anxF = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.gender = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.interaction = calc_HDI((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.intercept = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.estimates.anxF = FDT.fit.SAG.mratio.SAGIS.modelFit.mu_beta1;
FDT.fit.SAG.mratio.SAGIS.estimates.anxM = FDT.fit.SAG.mratio.SAGIS.modelFit.mu_beta2;
FDT.fit.SAG.mratio.SAGIS.estimates.gender = FDT.fit.SAG.mratio.SAGIS.modelFit.mu_beta3;
FDT.fit.SAG.mratio.SAGIS.estimates.interaction = FDT.fit.SAG.mratio.SAGIS.modelFit.mu_beta1 - FDT.fit.SAG.mratio.SAGIS.modelFit.mu_beta2;

% Fit the linear model for sensitivity (filter number)
filters_SAG = FDT.data.summary.avgFilter;
filters_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG.sensitivity.SAGI.modelFit = fitlm(FDT.GLMs.SAGI, filters_SAG, 'CategoricalVars', 2);
FDT.fit.SAG.sensitivity.SAGI.coefficients = FDT.fit.SAG.sensitivity.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensitivity.SAGI.pvalues = FDT.fit.SAG.sensitivity.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG.sensitivity.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS, filters_SAG, 'CategoricalVars', 3);
FDT.fit.SAG.sensitivity.SAGIS.coefficients = FDT.fit.SAG.sensitivity.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensitivity.SAGIS.pvalues = FDT.fit.SAG.sensitivity.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_SAG = FDT.data.summary.avgConfidence;
confidence_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG.sensibility.SAGI.modelFit = fitlm(FDT.GLMs.SAGI, confidence_SAG, 'CategoricalVars', 2);
FDT.fit.SAG.sensibility.SAGI.coefficients = FDT.fit.SAG.sensibility.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensibility.SAGI.pvalues = FDT.fit.SAG.sensibility.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG.sensibility.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS, confidence_SAG, 'CategoricalVars', 3);
FDT.fit.SAG.sensibility.SAGIS.coefficients = FDT.fit.SAG.sensibility.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensibility.SAGIS.pvalues = FDT.fit.SAG.sensibility.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.SAG.bias.SAGI.modelFit = fitlm(FDT.GLMs.SAGI, FDT.fit.SAG.mratio.SAGI.modelFit.c1, 'CategoricalVars', 2);
FDT.fit.SAG.bias.SAGI.coefficients = FDT.fit.SAG.bias.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.bias.SAGI.pvalues = FDT.fit.SAG.bias.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG.bias.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS, FDT.fit.SAG.mratio.SAGIS.modelFit.c1, 'CategoricalVars', 3);
FDT.fit.SAG.bias.SAGIS.coefficients = FDT.fit.SAG.bias.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.bias.SAGIS.pvalues = FDT.fit.SAG.bias.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.SAG.dPrime.SAGI.modelFit = fitlm(FDT.GLMs.SAGI, FDT.fit.SAG.mratio.SAGI.modelFit.d1, 'CategoricalVars', 2);
FDT.fit.SAG.dPrime.SAGI.coefficients = FDT.fit.SAG.dPrime.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.dPrime.SAGI.pvalues = FDT.fit.SAG.dPrime.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG.dPrime.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS, FDT.fit.SAG.mratio.SAGI.modelFit.d1, 'CategoricalVars', 3);
FDT.fit.SAG.dPrime.SAGIS.coefficients = FDT.fit.SAG.dPrime.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.dPrime.SAGIS.pvalues = FDT.fit.SAG.dPrime.SAGIS.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE REDUCED SA/SAG MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create reduced SA GLM
FDT.GLMs.SA_red = FDT.GLMs.SA(:,1);

% Fit the HMeta-d model for mratio
FDT.fit.SA_red.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SA.nR_S1, FDT.data.trials2counts.SA.nR_S2, FDT.GLMs.SA_red', FDT.settings.fit.mcmc_params);
FDT.fit.SA_red.mratio.HDI.anx_corr = calc_HDI(FDT.fit.SA_red.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.SA_red.mratio.HDI.anx = calc_HDI(FDT.fit.SA_red.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SA_red.mratio.HDI.intercept = calc_HDI(FDT.fit.SA_red.mratio.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
% Plot the HDI for checking
figure;
set(gcf, 'Position', [200 200 400 300])
histogram(FDT.fit.SA_red.mratio.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.SA_red.mratio.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SA_red.mratio.HDI.anx(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SA_red.mratio.HDI.anx(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('State anxiety Beta HDI');
set(gca, 'FontSize', 14)

% Fit the linear model for sensitivity (filter number)
filters_SA = FDT.data.summary.avgFilter;
filters_SA(idxNaN_row_SA) = [];
FDT.fit.SA_red.sensitivity.modelFit = fitlm(FDT.GLMs.SA_red, filters_SA);
FDT.fit.SA_red.sensitivity.coefficients = FDT.fit.SA_red.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.SA_red.sensitivity.pvalues = FDT.fit.SA_red.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_SA = FDT.data.summary.avgConfidence;
confidence_SA(idxNaN_row_SA) = [];
FDT.fit.SA_red.sensibility.modelFit = fitlm(FDT.GLMs.SA_red, confidence_SA);
FDT.fit.SA_red.sensibility.coefficients = FDT.fit.SA_red.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.SA_red.sensibility.pvalues = FDT.fit.SA_red.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.SA_red.bias.modelFit = fitlm(FDT.GLMs.SA_red, FDT.fit.SA_red.mratio.modelFit.c1);
FDT.fit.SA_red.bias.coefficients = FDT.fit.SA_red.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.SA_red.bias.pvalues = FDT.fit.SA_red.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.SA_red.dPrime.modelFit = fitlm(FDT.GLMs.SA_red, FDT.fit.SA_red.mratio.modelFit.d1);
FDT.fit.SA_red.dPrime.coefficients = FDT.fit.SA_red.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.SA_red.dPrime.pvalues = FDT.fit.SA_red.dPrime.modelFit.Coefficients.pValue(:);

% Create reduced SAG GLMs
FDT.GLMs.SAGI_red = FDT.GLMs.SAGI(:,1:3);
FDT.GLMs.SAGIS_red = FDT.GLMs.SAGIS(:,1:3);

% Fit the reduced HMeta-d model for mratio (specified interaction term)
FDT.fit.SAG_red.mratio.SAGI.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAG.nR_S1, FDT.data.trials2counts.SAG.nR_S2, FDT.GLMs.SAGI_red', FDT.settings.fit.mcmc_params);
FDT.fit.SAG_red.mratio.SAGI.HDI.anxF = calc_HDI((FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:)), 0.95);
FDT.fit.SAG_red.mratio.SAGI.HDI.anxM = calc_HDI(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SAG_red.mratio.SAGI.HDI.gender = calc_HDI(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG_red.mratio.SAGI.HDI.interaction = calc_HDI(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.SAG_red.mratio.SAGI.HDI.intercept = calc_HDI(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.SAG_red.mratio.SAGI.estimates.anxF = FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta3;
FDT.fit.SAG_red.mratio.SAGI.estimates.anxM = FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta1;
FDT.fit.SAG_red.mratio.SAGI.estimates.gender = FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta2;
FDT.fit.SAG_red.mratio.SAGI.estimates.interaction = FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta3;
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 800 300])
subplot(1,2,1)
histogram(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_logMratio(:));
hold on
xline(FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_logMratio, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.intercept(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.intercept(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Intercept (Male) HDI');
set(gca, 'FontSize', 14)
subplot(1,2,2)
histogram(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:));
hold on
xline(FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta2, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.gender(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.gender(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Gender (F>M) Beta HDI');
set(gca, 'FontSize', 14)
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 1200 300])
subplot(1,3,1)
histogram((FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:)));
hold on
xline((FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta1 + FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta3), 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.anxF(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.anxF(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Female State anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,2)
histogram(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.anxM(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.anxM(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Male State anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,3)
histogram(FDT.fit.SAG_red.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:));
hold on
xline(FDT.fit.SAG_red.mratio.SAGI.modelFit.mu_beta3, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.interaction(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.SAG_red.mratio.SAGI.HDI.interaction(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Interaction Beta HDI');
set(gca, 'FontSize', 14)

% Fit the reduced HMeta-d model for mratio (split interaction term)
FDT.fit.SAG_red.mratio.SAGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAG.nR_S1, FDT.data.trials2counts.SAG.nR_S2, FDT.GLMs.SAGIS_red', FDT.settings.fit.mcmc_params);
FDT.fit.SAG_red.mratio.SAGIS.HDI.anxF = calc_HDI(FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SAG_red.mratio.SAGIS.HDI.anxM = calc_HDI(FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG_red.mratio.SAGIS.HDI.gender = calc_HDI(FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.SAG_red.mratio.SAGIS.HDI.interaction = calc_HDI((FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.SAG_red.mratio.SAGIS.HDI.intercept = calc_HDI(FDT.fit.SAG_red.mratio.SAGIS.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.SAG_red.mratio.SAGIS.estimates.anxF = FDT.fit.SAG_red.mratio.SAGIS.modelFit.mu_beta1;
FDT.fit.SAG_red.mratio.SAGIS.estimates.anxM = FDT.fit.SAG_red.mratio.SAGIS.modelFit.mu_beta2;
FDT.fit.SAG_red.mratio.SAGIS.estimates.gender = FDT.fit.SAG_red.mratio.SAGIS.modelFit.mu_beta3;
FDT.fit.SAG_red.mratio.SAGIS.estimates.interaction = FDT.fit.SAG_red.mratio.SAGIS.modelFit.mu_beta1 - FDT.fit.SAG_red.mratio.SAGIS.modelFit.mu_beta2;

% Fit the linear model for sensitivity (filter number)
filters_SAG = FDT.data.summary.avgFilter;
filters_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG_red.sensitivity.SAGI.modelFit = fitlm(FDT.GLMs.SAGI_red, filters_SAG, 'CategoricalVars', 2);
FDT.fit.SAG_red.sensitivity.SAGI.coefficients = FDT.fit.SAG_red.sensitivity.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.sensitivity.SAGI.pvalues = FDT.fit.SAG_red.sensitivity.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG_red.sensitivity.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS_red, filters_SAG, 'CategoricalVars', 3);
FDT.fit.SAG_red.sensitivity.SAGIS.coefficients = FDT.fit.SAG_red.sensitivity.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.sensitivity.SAGIS.pvalues = FDT.fit.SAG_red.sensitivity.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_SAG = FDT.data.summary.avgConfidence;
confidence_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG_red.sensibility.SAGI.modelFit = fitlm(FDT.GLMs.SAGI_red, confidence_SAG, 'CategoricalVars', 2);
FDT.fit.SAG_red.sensibility.SAGI.coefficients = FDT.fit.SAG_red.sensibility.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.sensibility.SAGI.pvalues = FDT.fit.SAG_red.sensibility.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG_red.sensibility.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS_red, confidence_SAG, 'CategoricalVars', 3);
FDT.fit.SAG_red.sensibility.SAGIS.coefficients = FDT.fit.SAG_red.sensibility.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.sensibility.SAGIS.pvalues = FDT.fit.SAG_red.sensibility.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.SAG_red.bias.SAGI.modelFit = fitlm(FDT.GLMs.SAGI_red, FDT.fit.SAG_red.mratio.SAGI.modelFit.c1, 'CategoricalVars', 2);
FDT.fit.SAG_red.bias.SAGI.coefficients = FDT.fit.SAG_red.bias.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.bias.SAGI.pvalues = FDT.fit.SAG_red.bias.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG_red.bias.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS_red, FDT.fit.SAG_red.mratio.SAGIS.modelFit.c1, 'CategoricalVars', 3);
FDT.fit.SAG_red.bias.SAGIS.coefficients = FDT.fit.SAG_red.bias.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.bias.SAGIS.pvalues = FDT.fit.SAG_red.bias.SAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.SAG_red.dPrime.SAGI.modelFit = fitlm(FDT.GLMs.SAGI_red, FDT.fit.SAG_red.mratio.SAGI.modelFit.d1, 'CategoricalVars', 2);
FDT.fit.SAG_red.dPrime.SAGI.coefficients = FDT.fit.SAG_red.dPrime.SAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.dPrime.SAGI.pvalues = FDT.fit.SAG_red.dPrime.SAGI.modelFit.Coefficients.pValue(:);
FDT.fit.SAG_red.dPrime.SAGIS.modelFit = fitlm(FDT.GLMs.SAGIS_red, FDT.fit.SAG_red.mratio.SAGI.modelFit.d1, 'CategoricalVars', 3);
FDT.fit.SAG_red.dPrime.SAGIS.coefficients = FDT.fit.SAG_red.dPrime.SAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG_red.dPrime.SAGIS.pvalues = FDT.fit.SAG_red.dPrime.SAGIS.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE D GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the depression variable
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));

% Create the GLM matrix
Dep = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
end
FDT.GLMs.D(:,1) = Dep;
FDT.GLMs.D(1:length(FDT.data.IDs{1}),2:4) = -1;
FDT.GLMs.D(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),2) = 1;
FDT.GLMs.D(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),3) = 1;
FDT.GLMs.D(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),4) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_D col] = find(isnan(FDT.GLMs.D));
FDT.GLMs.D(idxNaN_row_D,:) = [];
FDT.data.trials2counts.D.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.D.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.D.nR_S1(idxNaN_row_D) = [];
FDT.data.trials2counts.D.nR_S2(idxNaN_row_D) = [];

% Zscore variables
FDT.GLMs.D(:,1) = zscore(FDT.GLMs.D(:,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE D MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio
FDT.fit.D.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.D.nR_S1, FDT.data.trials2counts.D.nR_S2, FDT.GLMs.D', FDT.settings.fit.mcmc_params);
FDT.fit.D.mratio.HDI.dep_corr = calc_HDI(FDT.fit.D.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.D.mratio.HDI.dep = calc_HDI(FDT.fit.D.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.D.mratio.HDI.intercept = calc_HDI(FDT.fit.D.mratio.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
% Plot the HDI for checking
figure;
set(gcf, 'Position', [200 200 400 300])
histogram(FDT.fit.D.mratio.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.D.mratio.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.D.mratio.HDI.dep(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.D.mratio.HDI.dep(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Depression Beta HDI');
set(gca, 'FontSize', 14)

% Fit the linear model for sensitivity (filter number)
filters_D = FDT.data.summary.avgFilter;
filters_D(idxNaN_row_D) = [];
FDT.fit.D.sensitivity.modelFit = fitlm(FDT.GLMs.D, filters_D);
FDT.fit.D.sensitivity.coefficients = FDT.fit.D.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.D.sensitivity.pvalues = FDT.fit.D.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_D = FDT.data.summary.avgConfidence;
confidence_D(idxNaN_row_D) = [];
FDT.fit.D.sensibility.modelFit = fitlm(FDT.GLMs.D, confidence_D);
FDT.fit.D.sensibility.coefficients = FDT.fit.D.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.D.sensibility.pvalues = FDT.fit.D.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.D.bias.modelFit = fitlm(FDT.GLMs.D, FDT.fit.D.mratio.modelFit.c1);
FDT.fit.D.bias.coefficients = FDT.fit.D.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.D.bias.pvalues = FDT.fit.D.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.D.dPrime.modelFit = fitlm(FDT.GLMs.D, FDT.fit.D.mratio.modelFit.d1);
FDT.fit.D.dPrime.coefficients = FDT.fit.D.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.D.dPrime.pvalues = FDT.fit.D.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE DG GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the depression variable
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));

% Find the gender variable
idxGen = find(strcmp(FDT.settings.variables.summarySheet, 'GENDER'));

% Create the gender variable and GLM with interaction specified in terms matrix
Dep = [];
Gender = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
    Gender = [Gender; FDT.data.raw.extra{study}(:,idxGen)];
end
FDT.GLMs.DG(:,1) = Dep;
FDT.GLMs.DG(:,2) = Gender;
FDT.GLMs.DG(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.DG(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.DG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.DG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_DG col] = find(isnan(FDT.GLMs.DG));
FDT.GLMs.DG(idxNaN_row_DG,:) = [];
FDT.data.trials2counts.DG.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.DG.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.DG.nR_S1(idxNaN_row_DG) = [];
FDT.data.trials2counts.DG.nR_S2(idxNaN_row_DG) = [];

% Zscore variables
FDT.GLMs.DG(:,1) = zscore(FDT.GLMs.DG(:,1));

% Create the SAG matrix with an explicit interaction term
FDT.GLMs.DGI = FDT.GLMs.DG;
FDT.GLMs.DGI(:,4:6) = FDT.GLMs.DG(:,3:5);
FDT.GLMs.DGI(:,3) = FDT.GLMs.DG(:,1) .* FDT.GLMs.DG(:,2);

% Create the SAG matrix with a split interaction term
FDT.GLMs.DGIS(:,1) = FDT.GLMs.DG(:,1);
FDT.GLMs.DGIS(:,2) = FDT.GLMs.DG(:,1);
FDT.GLMs.DGIS(:,3) = FDT.GLMs.DG(:,2);
FDT.GLMs.DGIS((FDT.GLMs.DG(:,2) == 0),1) = 0;
FDT.GLMs.DGIS((FDT.GLMs.DG(:,2) == 1),2) = 0;
FDT.GLMs.DGIS(:,4:6) = FDT.GLMs.DG(:,3:5);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE DG MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio (specified interaction term)
FDT.fit.DG.mratio.DGI.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.DG.nR_S1, FDT.data.trials2counts.DG.nR_S2, FDT.GLMs.DGI', FDT.settings.fit.mcmc_params);
FDT.fit.DG.mratio.DGI.HDI.depF = calc_HDI((FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta3(:)), 0.95);
FDT.fit.DG.mratio.DGI.HDI.depM = calc_HDI(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.DG.mratio.DGI.HDI.gender = calc_HDI(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.DG.mratio.DGI.HDI.interaction = calc_HDI(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.DG.mratio.DGI.HDI.intercept = calc_HDI(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.DG.mratio.DGI.estimates.depF = FDT.fit.DG.mratio.DGI.modelFit.mu_beta1 + FDT.fit.DG.mratio.DGI.modelFit.mu_beta3;
FDT.fit.DG.mratio.DGI.estimates.depM = FDT.fit.DG.mratio.DGI.modelFit.mu_beta1;
FDT.fit.DG.mratio.DGI.estimates.gender = FDT.fit.DG.mratio.DGI.modelFit.mu_beta2;
FDT.fit.DG.mratio.DGI.estimates.interaction = FDT.fit.DG.mratio.DGI.modelFit.mu_beta3;
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 800 300])
subplot(1,2,1)
histogram(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_logMratio(:));
hold on
xline(FDT.fit.DG.mratio.DGI.modelFit.mu_logMratio, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.intercept(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.intercept(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Intercept (Male) HDI');
set(gca, 'FontSize', 14)
subplot(1,2,2)
histogram(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta2(:));
hold on
xline(FDT.fit.DG.mratio.DGI.modelFit.mu_beta2, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.gender(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.gender(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Gender (F>M) Beta HDI');
set(gca, 'FontSize', 14)
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 1200 300])
subplot(1,3,1)
histogram((FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta3(:)));
hold on
xline((FDT.fit.DG.mratio.DGI.modelFit.mu_beta1 + FDT.fit.DG.mratio.DGI.modelFit.mu_beta3), 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.depF(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.depF(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Female Depression Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,2)
histogram(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.DG.mratio.DGI.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.depM(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.depM(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Male Depression Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,3)
histogram(FDT.fit.DG.mratio.DGI.modelFit.mcmc.samples.mu_beta3(:));
hold on
xline(FDT.fit.DG.mratio.DGI.modelFit.mu_beta3, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.interaction(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.DG.mratio.DGI.HDI.interaction(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Interaction Beta HDI');
set(gca, 'FontSize', 14)

% Fit the HMeta-d model for mratio (split interaction term)
FDT.fit.DG.mratio.DGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.DG.nR_S1, FDT.data.trials2counts.DG.nR_S2, FDT.GLMs.DGIS', FDT.settings.fit.mcmc_params);
FDT.fit.DG.mratio.DGIS.HDI.depF = calc_HDI(FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.DG.mratio.DGIS.HDI.depM = calc_HDI(FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.DG.mratio.DGIS.HDI.gender = calc_HDI(FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.DG.mratio.DGIS.HDI.interaction = calc_HDI((FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.DG.mratio.DGIS.HDI.intercept = calc_HDI(FDT.fit.DG.mratio.DGIS.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.DG.mratio.DGIS.estimates.depF = FDT.fit.DG.mratio.DGIS.modelFit.mu_beta1;
FDT.fit.DG.mratio.DGIS.estimates.depM = FDT.fit.DG.mratio.DGIS.modelFit.mu_beta2;
FDT.fit.DG.mratio.DGIS.estimates.gender = FDT.fit.DG.mratio.DGIS.modelFit.mu_beta3;
FDT.fit.DG.mratio.DGIS.estimates.interaction = FDT.fit.DG.mratio.DGIS.modelFit.mu_beta1 - FDT.fit.DG.mratio.DGIS.modelFit.mu_beta2;

% Fit the linear model for sensitivity (filter number)
filters_DG = FDT.data.summary.avgFilter;
filters_DG(idxNaN_row_DG) = [];
FDT.fit.DG.sensitivity.DGI.modelFit = fitlm(FDT.GLMs.DGI, filters_DG, 'CategoricalVars', 2);
FDT.fit.DG.sensitivity.DGI.coefficients = FDT.fit.DG.sensitivity.DGI.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.sensitivity.DGI.pvalues = FDT.fit.DG.sensitivity.DGI.modelFit.Coefficients.pValue(:);
FDT.fit.DG.sensitivity.DGIS.modelFit = fitlm(FDT.GLMs.DGIS, filters_DG, 'CategoricalVars', 3);
FDT.fit.DG.sensitivity.DGIS.coefficients = FDT.fit.DG.sensitivity.DGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.sensitivity.DGIS.pvalues = FDT.fit.DG.sensitivity.DGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_DG = FDT.data.summary.avgConfidence;
confidence_DG(idxNaN_row_DG) = [];
FDT.fit.DG.sensibility.DGI.modelFit = fitlm(FDT.GLMs.DGI, confidence_DG, 'CategoricalVars', 2);
FDT.fit.DG.sensibility.DGI.coefficients = FDT.fit.DG.sensibility.DGI.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.sensibility.DGI.pvalues = FDT.fit.DG.sensibility.DGI.modelFit.Coefficients.pValue(:);
FDT.fit.DG.sensibility.DGIS.modelFit = fitlm(FDT.GLMs.DGIS, confidence_DG, 'CategoricalVars', 3);
FDT.fit.DG.sensibility.DGIS.coefficients = FDT.fit.DG.sensibility.DGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.sensibility.DGIS.pvalues = FDT.fit.DG.sensibility.DGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.DG.bias.DGI.modelFit = fitlm(FDT.GLMs.DGI, FDT.fit.DG.mratio.DGI.modelFit.c1, 'CategoricalVars', 2);
FDT.fit.DG.bias.DGI.coefficients = FDT.fit.DG.bias.DGI.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.bias.DGI.pvalues = FDT.fit.DG.bias.DGI.modelFit.Coefficients.pValue(:);
FDT.fit.DG.bias.DGIS.modelFit = fitlm(FDT.GLMs.DGIS, FDT.fit.DG.mratio.DGIS.modelFit.c1, 'CategoricalVars', 3);
FDT.fit.DG.bias.DGIS.coefficients = FDT.fit.DG.bias.DGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.bias.DGIS.pvalues = FDT.fit.DG.bias.DGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.DG.dPrime.DGI.modelFit = fitlm(FDT.GLMs.DGI, FDT.fit.DG.mratio.DGI.modelFit.d1, 'CategoricalVars', 2);
FDT.fit.DG.dPrime.DGI.coefficients = FDT.fit.DG.dPrime.DGI.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.dPrime.DGI.pvalues = FDT.fit.DG.dPrime.DGI.modelFit.Coefficients.pValue(:);
FDT.fit.DG.dPrime.DGIS.modelFit = fitlm(FDT.GLMs.DGIS, FDT.fit.DG.mratio.DGI.modelFit.d1, 'CategoricalVars', 3);
FDT.fit.DG.dPrime.DGIS.coefficients = FDT.fit.DG.dPrime.DGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.DG.dPrime.DGIS.pvalues = FDT.fit.DG.dPrime.DGIS.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE TA GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the trait anxiety variable
idxTA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-T'));

% Create the GLM for state anxiety (with regressors for study location)
TA = [];
for study = 1:4
    TA = [TA; FDT.data.raw.extra{study}(:,idxTA)];
end
FDT.GLMs.TA = [TA, zeros(length(TA),3)];
FDT.GLMs.TA(1:length(FDT.data.IDs{1}),2:4) = -1;
FDT.GLMs.TA(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),2) = 1;
FDT.GLMs.TA(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),3) = 1;
FDT.GLMs.TA(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),4) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_TA col] = find(isnan(TA));
FDT.GLMs.TA(idxNaN_row_TA,:) = [];
FDT.data.trials2counts.TA.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.TA.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.TA.nR_S1(idxNaN_row_TA) = [];
FDT.data.trials2counts.TA.nR_S2(idxNaN_row_TA) = [];

% Zscore variables
FDT.GLMs.TA(:,1) = zscore(FDT.GLMs.TA(:,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE TA MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio
FDT.fit.TA.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TA.nR_S1, FDT.data.trials2counts.TA.nR_S2, FDT.GLMs.TA', FDT.settings.fit.mcmc_params);
FDT.fit.TA.mratio.HDI.anx_corr = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.TA.mratio.HDI.anx = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.TA.mratio.HDI.intercept = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
% Plot the HDI for checking
figure;
set(gcf, 'Position', [200 200 400 300])
histogram(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.TA.mratio.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TA.mratio.HDI.anx(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TA.mratio.HDI.anx(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Trait anxiety Beta HDI');
set(gca, 'FontSize', 14)

% Fit the linear model for sensitivity (filter number)
filters_TA = FDT.data.summary.avgFilter;
filters_TA(idxNaN_row_TA) = [];
FDT.fit.TA.sensitivity.modelFit = fitlm(FDT.GLMs.TA, filters_TA);
FDT.fit.TA.sensitivity.coefficients = FDT.fit.TA.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.TA.sensitivity.pvalues = FDT.fit.TA.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_TA = FDT.data.summary.avgConfidence;
confidence_TA(idxNaN_row_TA) = [];
FDT.fit.TA.sensibility.modelFit = fitlm(FDT.GLMs.TA, confidence_TA);
FDT.fit.TA.sensibility.coefficients = FDT.fit.TA.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.TA.sensibility.pvalues = FDT.fit.TA.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.TA.bias.modelFit = fitlm(FDT.GLMs.TA, FDT.fit.TA.mratio.modelFit.c1);
FDT.fit.TA.bias.coefficients = FDT.fit.TA.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.TA.bias.pvalues = FDT.fit.TA.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.TA.dPrime.modelFit = fitlm(FDT.GLMs.TA, FDT.fit.TA.mratio.modelFit.d1);
FDT.fit.TA.dPrime.coefficients = FDT.fit.TA.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.TA.dPrime.pvalues = FDT.fit.TA.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE TAG GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the trait anxiety variable
idxTA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-T'));

% Find the gender variable
idxGen = find(strcmp(FDT.settings.variables.summarySheet, 'GENDER'));

% Create the gender variable and GLM with interaction specified in terms matrix
TA = [];
Gender = [];
for study = 1:4
    TA = [TA; FDT.data.raw.extra{study}(:,idxTA)];
    Gender = [Gender; FDT.data.raw.extra{study}(:,idxGen)];
end
FDT.GLMs.TAG(:,1) = TA;
FDT.GLMs.TAG(:,2) = Gender;
FDT.GLMs.TAG(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_TAG col] = find(isnan(FDT.GLMs.TAG));
FDT.GLMs.TAG(idxNaN_row_TAG,:) = [];
FDT.data.trials2counts.TAG.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.TAG.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.TAG.nR_S1(idxNaN_row_TAG) = [];
FDT.data.trials2counts.TAG.nR_S2(idxNaN_row_TAG) = [];

% Zscore variables
FDT.GLMs.TAG(:,1) = zscore(FDT.GLMs.TAG(:,1));

% Create the TAG matrix with an explicit interaction term
FDT.GLMs.TAGI = FDT.GLMs.TAG;
FDT.GLMs.TAGI(:,4:6) = FDT.GLMs.TAG(:,3:5);
FDT.GLMs.TAGI(:,3) = FDT.GLMs.TAG(:,1) .* FDT.GLMs.TAG(:,2);

% Create the TAG matrix with a split interaction term
FDT.GLMs.TAGIS(:,1) = FDT.GLMs.TAG(:,1);
FDT.GLMs.TAGIS(:,2) = FDT.GLMs.TAG(:,1);
FDT.GLMs.TAGIS(:,3) = FDT.GLMs.TAG(:,2);
FDT.GLMs.TAGIS((FDT.GLMs.TAG(:,2) == 0),1) = 0;
FDT.GLMs.TAGIS((FDT.GLMs.TAG(:,2) == 1),2) = 0;
FDT.GLMs.TAGIS(:,4:6) = FDT.GLMs.TAG(:,3:5);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE TAG MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio (specified interaction term)
FDT.fit.TAG.mratio.TAGI.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TAG.nR_S1, FDT.data.trials2counts.TAG.nR_S2, FDT.GLMs.TAGI', FDT.settings.fit.mcmc_params);
FDT.fit.TAG.mratio.TAGI.HDI.anxF = calc_HDI((FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:)), 0.95);
FDT.fit.TAG.mratio.TAGI.HDI.anxM = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.TAG.mratio.TAGI.HDI.gender = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.TAG.mratio.TAGI.HDI.interaction = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.TAG.mratio.TAGI.HDI.intercept = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.TAG.mratio.TAGI.estimates.anxF = FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1 + FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3;
FDT.fit.TAG.mratio.TAGI.estimates.anxM = FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1;
FDT.fit.TAG.mratio.TAGI.estimates.gender = FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta2;
FDT.fit.TAG.mratio.TAGI.estimates.interaction = FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3;
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 800 300])
subplot(1,2,1)
histogram(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_logMratio(:));
hold on
xline(FDT.fit.TAG.mratio.TAGI.modelFit.mu_logMratio, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.intercept(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.intercept(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Intercept (Male) HDI');
set(gca, 'FontSize', 14)
subplot(1,2,2)
histogram(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:));
hold on
xline(FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta2, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.gender(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.gender(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Gender (F>M) Beta HDI');
set(gca, 'FontSize', 14)
% Plot the gender HDIs for checking
figure;
set(gcf, 'Position', [200 200 1200 300])
subplot(1,3,1)
histogram((FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:)));
hold on
xline((FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1 + FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3), 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.anxF(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.anxF(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Female Trait anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,2)
histogram(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:));
hold on
xline(FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta1, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.anxM(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.anxM(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Male Trait anxiety Beta HDI');
set(gca, 'FontSize', 14)
subplot(1,3,3)
histogram(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:));
hold on
xline(FDT.fit.TAG.mratio.TAGI.modelFit.mu_beta3, 'color', 'b', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.interaction(1), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(FDT.fit.TAG.mratio.TAGI.HDI.interaction(2), 'linestyle', '--', 'color', 'k', 'linewidth', 1.5);
xline(0, 'linestyle', ':', 'color', 'r', 'linewidth', 2);
xlabel('Parameter');
ylabel('Sample count');
title('Interaction Beta HDI');
set(gca, 'FontSize', 14)

% Fit the HMeta-d model for mratio (split interaction term)
FDT.fit.TAG.mratio.TAGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TAG.nR_S1, FDT.data.trials2counts.TAG.nR_S2, FDT.GLMs.TAGIS', FDT.settings.fit.mcmc_params);
FDT.fit.TAG.mratio.TAGIS.HDI.anxF = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.TAG.mratio.TAGIS.HDI.anxM = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.TAG.mratio.TAGIS.HDI.gender = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.TAG.mratio.TAGIS.HDI.interaction = calc_HDI((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.TAG.mratio.TAGIS.HDI.intercept = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_logMratio(:), 0.95);
FDT.fit.TAG.mratio.TAGIS.estimates.anxF = FDT.fit.TAG.mratio.TAGIS.modelFit.mu_beta1;
FDT.fit.TAG.mratio.TAGIS.estimates.anxM = FDT.fit.TAG.mratio.TAGIS.modelFit.mu_beta2;
FDT.fit.TAG.mratio.TAGIS.estimates.gender = FDT.fit.TAG.mratio.TAGIS.modelFit.mu_beta3;
FDT.fit.TAG.mratio.TAGIS.estimates.interaction = FDT.fit.TAG.mratio.TAGIS.modelFit.mu_beta1 - FDT.fit.TAG.mratio.TAGIS.modelFit.mu_beta2;

% Fit the linear model for sensitivity (filter number)
filters_TAG = FDT.data.summary.avgFilter;
filters_TAG(idxNaN_row_TAG) = [];
FDT.fit.TAG.sensitivity.TAGI.modelFit = fitlm(FDT.GLMs.TAGI, filters_TAG, 'CategoricalVars', 2);
FDT.fit.TAG.sensitivity.TAGI.coefficients = FDT.fit.TAG.sensitivity.TAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.sensitivity.TAGI.pvalues = FDT.fit.TAG.sensitivity.TAGI.modelFit.Coefficients.pValue(:);
FDT.fit.TAG.sensitivity.TAGIS.modelFit = fitlm(FDT.GLMs.TAGIS, filters_TAG, 'CategoricalVars', 3);
FDT.fit.TAG.sensitivity.TAGIS.coefficients = FDT.fit.TAG.sensitivity.TAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.sensitivity.TAGIS.pvalues = FDT.fit.TAG.sensitivity.TAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_TAG = FDT.data.summary.avgConfidence;
confidence_TAG(idxNaN_row_TAG) = [];
FDT.fit.TAG.sensibility.TAGI.modelFit = fitlm(FDT.GLMs.TAGI, confidence_TAG, 'CategoricalVars', 2);
FDT.fit.TAG.sensibility.TAGI.coefficients = FDT.fit.TAG.sensibility.TAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.sensibility.TAGI.pvalues = FDT.fit.TAG.sensibility.TAGI.modelFit.Coefficients.pValue(:);
FDT.fit.TAG.sensibility.TAGIS.modelFit = fitlm(FDT.GLMs.TAGIS, confidence_TAG, 'CategoricalVars', 3);
FDT.fit.TAG.sensibility.TAGIS.coefficients = FDT.fit.TAG.sensibility.TAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.sensibility.TAGIS.pvalues = FDT.fit.TAG.sensibility.TAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.TAG.bias.TAGI.modelFit = fitlm(FDT.GLMs.TAGI, FDT.fit.TAG.mratio.TAGI.modelFit.c1, 'CategoricalVars', 2);
FDT.fit.TAG.bias.TAGI.coefficients = FDT.fit.TAG.bias.TAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.bias.TAGI.pvalues = FDT.fit.TAG.bias.TAGI.modelFit.Coefficients.pValue(:);
FDT.fit.TAG.bias.TAGIS.modelFit = fitlm(FDT.GLMs.TAGIS, FDT.fit.TAG.mratio.TAGIS.modelFit.c1, 'CategoricalVars', 3);
FDT.fit.TAG.bias.TAGIS.coefficients = FDT.fit.TAG.bias.TAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.bias.TAGIS.pvalues = FDT.fit.TAG.bias.TAGIS.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.TAG.dPrime.TAGI.modelFit = fitlm(FDT.GLMs.TAGI, FDT.fit.TAG.mratio.TAGI.modelFit.d1, 'CategoricalVars', 2);
FDT.fit.TAG.dPrime.TAGI.coefficients = FDT.fit.TAG.dPrime.TAGI.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.dPrime.TAGI.pvalues = FDT.fit.TAG.dPrime.TAGI.modelFit.Coefficients.pValue(:);
FDT.fit.TAG.dPrime.TAGIS.modelFit = fitlm(FDT.GLMs.TAGIS, FDT.fit.TAG.mratio.TAGI.modelFit.d1, 'CategoricalVars', 3);
FDT.fit.TAG.dPrime.TAGIS.coefficients = FDT.fit.TAG.dPrime.TAGIS.modelFit.Coefficients.Estimate(:);
FDT.fit.TAG.dPrime.TAGIS.pvalues = FDT.fit.TAG.dPrime.TAGIS.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE OVERALL SUMMARY MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accuracy
FDT.values.accuracy.mean = mean(FDT.data.summary.accuracy);
FDT.values.accuracy.sd = std(FDT.data.summary.accuracy);
FDT.values.accuracy.se = FDT.values.accuracy.sd / sqrt(length(FDT.data.summary.accuracy));

% Task difficulty (d')
FDT.values.dPrime.mean = mean(FDT.fit.SA.mratio.modelFit.d1);
FDT.values.dPrime.sd = std(FDT.fit.SA.mratio.modelFit.d1);
FDT.values.dPrime.se = FDT.values.dPrime.sd / sqrt(length(FDT.fit.SA.mratio.modelFit.d1));

% Number of trials
FDT.values.trials.mean = mean(FDT.data.summary.trials);
FDT.values.trials.sd = std(FDT.data.summary.trials);
FDT.values.trials.se = FDT.values.trials.sd / sqrt(length(FDT.data.summary.trials));

% Sensitivity (filter number)
FDT.values.sensitivity.mean = mean(FDT.data.summary.avgFilter);
FDT.values.sensitivity.sd = std(FDT.data.summary.avgFilter);
FDT.values.sensitivity.se = FDT.values.sensitivity.sd / sqrt(length(FDT.data.summary.avgFilter));

% Bias
FDT.values.bias.mean = mean(FDT.fit.SA.mratio.modelFit.c1);
FDT.values.bias.sd = std(FDT.fit.SA.mratio.modelFit.c1);
FDT.values.bias.se = FDT.values.bias.sd / sqrt(length(FDT.fit.SA.mratio.modelFit.c1));

% Sensibility (average confidence)
FDT.values.sensibility.mean = mean(FDT.data.summary.avgConfidence);
FDT.values.sensibility.sd = std(FDT.data.summary.avgConfidence);
FDT.values.sensibility.se = FDT.values.sensibility.sd / sqrt(length(FDT.data.summary.avgConfidence));

% logMratio
FDT.values.logMratio.mean = FDT.fit.SA.mratio.modelFit.mu_logMratio;
FDT.values.logMratio.sd = FDT.fit.SA.mratio.modelFit.sigma_logMratio;
FDT.values.logMratio.se = std(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_logMratio(:));

% Mratio
FDT.values.mratio.mean = exp(FDT.fit.SA.mratio.modelFit.mu_logMratio);
FDT.values.mratio.se = std(exp(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_logMratio(:)));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Remove samples from model fits to save space
% FDT.fit.SA.mratio.modelFit.mcmc.samples = [];
% FDT.fit.TA.mratio.modelFit.mcmc.samples = [];
% FDT.fit.TAD.mratio.modelFit.mcmc.samples = [];
% FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples = [];
% FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples = [];
% FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples = [];
% FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples = [];

% Save data matrix
save(FDT.settings.names.output, 'FDT', '-v7.3');


end

