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
FDT.fit.SA.mratio.HDIinterval = [0.0 0.987];
FDT.fit.SA.mratio.HDIinterval_unc = [0.0 0.95];
FDT.fit.SA.mratio.tails = 1;
FDT.fit.SA.mratio.HDI.anx_quantile = quantile(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SA.mratio.HDIinterval);
FDT.fit.SA.mratio.HDI.anx_inner = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
FDT.fit.SA.mratio.HDI.anx_outer = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SA.mratio.HDI.anx = [FDT.fit.SA.mratio.HDI.anx_outer(1) FDT.fit.SA.mratio.HDI.anx_inner(2)];
FDT.fit.SA.mratio.HDI.anx_quantile_unc = quantile(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SA.mratio.HDIinterval_unc);
FDT.fit.SA.mratio.HDI.anx_inner_unc = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9);
FDT.fit.SA.mratio.HDI.anx_outer_unc = calc_HDI(FDT.fit.SA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SA.mratio.HDI.anx_unc = [FDT.fit.SA.mratio.HDI.anx_outer_unc(1) FDT.fit.SA.mratio.HDI.anx_inner_unc(2)];

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
% CREATE THE TA GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the state anxiety variable
idxTA = find(strcmp(FDT.settings.variables.summarySheet, 'ANXIETY-T'));

% Create the GLM for trait anxiety (with regressors for study location)
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

% % Fit the HMeta-d model for mratio
% FDT.fit.TA.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TA.nR_S1, FDT.data.trials2counts.TA.nR_S2, FDT.GLMs.TA', FDT.settings.fit.mcmc_params);
% FDT.fit.TA.mratio.HDIinterval = [0.0 0.987];
% FDT.fit.TA.mratio.HDIinterval_unc = [0.0 0.95];
% FDT.fit.TA.mratio.tails = 1;
% FDT.fit.TA.mratio.HDI.anx_quantile = quantile(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TA.mratio.HDIinterval);
% FDT.fit.TA.mratio.HDI.anx_inner = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
% FDT.fit.TA.mratio.HDI.anx_outer = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.TA.mratio.HDI.anx = [FDT.fit.TA.mratio.HDI.anx_outer(1) FDT.fit.TA.mratio.HDI.anx_inner(2)];
% FDT.fit.TA.mratio.HDI.anx_quantile_unc = quantile(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TA.mratio.HDIinterval_unc);
% FDT.fit.TA.mratio.HDI.anx_inner_unc = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9);
% FDT.fit.TA.mratio.HDI.anx_outer_unc = calc_HDI(FDT.fit.TA.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.TA.mratio.HDI.anx_unc = [FDT.fit.TA.mratio.HDI.anx_outer_unc(1) FDT.fit.TA.mratio.HDI.anx_inner_unc(2)];
% 
% % Fit the linear model for sensitivity (filter number)
% filters_TA = FDT.data.summary.avgFilter;
% filters_TA(idxNaN_row_TA) = [];
% FDT.fit.TA.sensitivity.modelFit = fitlm(FDT.GLMs.TA, filters_TA);
% FDT.fit.TA.sensitivity.coefficients = FDT.fit.TA.sensitivity.modelFit.Coefficients.Estimate(:);
% FDT.fit.TA.sensitivity.pvalues = FDT.fit.TA.sensitivity.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for sensibility (confidence)
% confidence_TA = FDT.data.summary.avgConfidence;
% confidence_TA(idxNaN_row_TA) = [];
% FDT.fit.TA.sensibility.modelFit = fitlm(FDT.GLMs.TA, confidence_TA);
% FDT.fit.TA.sensibility.coefficients = FDT.fit.TA.sensibility.modelFit.Coefficients.Estimate(:);
% FDT.fit.TA.sensibility.pvalues = FDT.fit.TA.sensibility.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for decision bias
% FDT.fit.TA.bias.modelFit = fitlm(FDT.GLMs.TA, FDT.fit.TA.mratio.modelFit.c1);
% FDT.fit.TA.bias.coefficients = FDT.fit.TA.bias.modelFit.Coefficients.Estimate(:);
% FDT.fit.TA.bias.pvalues = FDT.fit.TA.bias.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for d' (sanity check)
% FDT.fit.TA.dPrime.modelFit = fitlm(FDT.GLMs.TA, FDT.fit.TA.mratio.modelFit.d1);
% FDT.fit.TA.dPrime.coefficients = FDT.fit.TA.dPrime.modelFit.Coefficients.Estimate(:);
% FDT.fit.TA.dPrime.pvalues = FDT.fit.TA.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE DEPRESSION GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the depression variable
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));

% Create the GLM for depression (with regressors for study location)
Dep = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
end
FDT.GLMs.Dep = [Dep, zeros(length(Dep),3)];
FDT.GLMs.Dep(1:length(FDT.data.IDs{1}),2:4) = -1;
FDT.GLMs.Dep(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),2) = 1;
FDT.GLMs.Dep(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),3) = 1;
FDT.GLMs.Dep(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),4) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_Dep col] = find(isnan(Dep));
FDT.GLMs.Dep(idxNaN_row_Dep,:) = [];
FDT.data.trials2counts.Dep.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.Dep.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.Dep.nR_S1(idxNaN_row_Dep) = [];
FDT.data.trials2counts.Dep.nR_S2(idxNaN_row_Dep) = [];

% Zscore variables
FDT.GLMs.Dep(:,1) = zscore(FDT.GLMs.Dep(:,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE DEPRESSION MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Fit the HMeta-d model for mratio
% FDT.fit.Dep.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.Dep.nR_S1, FDT.data.trials2counts.Dep.nR_S2, FDT.GLMs.Dep', FDT.settings.fit.mcmc_params);
% FDT.fit.Dep.mratio.HDIinterval = [0.0 0.987];
% FDT.fit.Dep.mratio.HDIinterval_unc = [0.0 0.95];
% FDT.fit.Dep.mratio.tails = 1;
% FDT.fit.Dep.mratio.HDI.anx_quantile = quantile(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.Dep.mratio.HDIinterval);
% FDT.fit.Dep.mratio.HDI.anx_inner = calc_HDI(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
% FDT.fit.Dep.mratio.HDI.anx_outer = calc_HDI(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.Dep.mratio.HDI.anx = [FDT.fit.Dep.mratio.HDI.anx_outer(1) FDT.fit.Dep.mratio.HDI.anx_inner(2)];
% FDT.fit.Dep.mratio.HDI.anx_quantile_unc = quantile(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.Dep.mratio.HDIinterval_unc);
% FDT.fit.Dep.mratio.HDI.anx_inner_unc = calc_HDI(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9);
% FDT.fit.Dep.mratio.HDI.anx_outer_unc = calc_HDI(FDT.fit.Dep.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.Dep.mratio.HDI.anx_unc = [FDT.fit.Dep.mratio.HDI.anx_outer_unc(1) FDT.fit.Dep.mratio.HDI.anx_inner_unc(2)];
% 
% % Fit the linear model for sensitivity (filter number)
% filters_Dep = FDT.data.summary.avgFilter;
% filters_Dep(idxNaN_row_Dep) = [];
% FDT.fit.Dep.sensitivity.modelFit = fitlm(FDT.GLMs.Dep, filters_Dep);
% FDT.fit.Dep.sensitivity.coefficients = FDT.fit.Dep.sensitivity.modelFit.Coefficients.Estimate(:);
% FDT.fit.Dep.sensitivity.pvalues = FDT.fit.Dep.sensitivity.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for sensibility (confidence)
% confidence_Dep = FDT.data.summary.avgConfidence;
% confidence_Dep(idxNaN_row_Dep) = [];
% FDT.fit.Dep.sensibility.modelFit = fitlm(FDT.GLMs.Dep, confidence_Dep);
% FDT.fit.Dep.sensibility.coefficients = FDT.fit.Dep.sensibility.modelFit.Coefficients.Estimate(:);
% FDT.fit.Dep.sensibility.pvalues = FDT.fit.Dep.sensibility.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for decision bias
% FDT.fit.Dep.bias.modelFit = fitlm(FDT.GLMs.Dep, FDT.fit.Dep.mratio.modelFit.c1);
% FDT.fit.Dep.bias.coefficients = FDT.fit.Dep.bias.modelFit.Coefficients.Estimate(:);
% FDT.fit.Dep.bias.pvalues = FDT.fit.Dep.bias.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for d' (sanity check)
% FDT.fit.Dep.dPrime.modelFit = fitlm(FDT.GLMs.Dep, FDT.fit.Dep.mratio.modelFit.d1);
% FDT.fit.Dep.dPrime.coefficients = FDT.fit.Dep.dPrime.modelFit.Coefficients.Estimate(:);
% FDT.fit.Dep.dPrime.pvalues = FDT.fit.Dep.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE SAD GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the depression variable
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));

% Create the depression variable and create GLM matrix
Dep = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
end
FDT.GLMs.SAD(:,1) = SA;
FDT.GLMs.SAD(:,2) = Dep;
FDT.GLMs.SAD(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.SAD(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.SAD(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.SAD(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_SAD col] = find(isnan(FDT.GLMs.SAD));
FDT.GLMs.SAD(idxNaN_row_SAD,:) = [];
FDT.data.trials2counts.SAD.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.SAD.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.SAD.nR_S1(idxNaN_row_SAD) = [];
FDT.data.trials2counts.SAD.nR_S2(idxNaN_row_SAD) = [];

% Zscore variables
FDT.GLMs.SAD(:,1) = zscore(FDT.GLMs.SAD(:,1));
FDT.GLMs.SAD(:,2) = zscore(FDT.GLMs.SAD(:,2));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE SAD MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Fit the HMeta-d model for mratio
% FDT.fit.SAD.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAD.nR_S1, FDT.data.trials2counts.SAD.nR_S2, FDT.GLMs.SAD', FDT.settings.fit.mcmc_params);
% FDT.fit.SAD.mratio.HDIinterval = [0.0 0.987];
% FDT.fit.SAD.mratio.HDIinterval_unc = [0.0 0.95];
% FDT.fit.SAD.mratio.tails = 1;
% FDT.fit.SAD.mratio.HDI.anx_quantile = quantile(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAD.mratio.HDIinterval);
% FDT.fit.SAD.mratio.HDI.anx_inner = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
% FDT.fit.SAD.mratio.HDI.anx_outer = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.SAD.mratio.HDI.anx = [FDT.fit.SAD.mratio.HDI.anx_outer(1) FDT.fit.SAD.mratio.HDI.anx_inner(2)];
% FDT.fit.SAD.mratio.HDI.anx_quantile_unc = quantile(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAD.mratio.HDIinterval_unc);
% FDT.fit.SAD.mratio.HDI.anx_inner_unc = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9);
% FDT.fit.SAD.mratio.HDI.anx_outer_unc = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
% FDT.fit.SAD.mratio.HDI.anx_unc = [FDT.fit.SAD.mratio.HDI.anx_outer_unc(1) FDT.fit.SAD.mratio.HDI.anx_inner_unc(2)];
% FDT.fit.SAD.mratio.HDI.dep_quantile = quantile(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAD.mratio.HDIinterval);
% FDT.fit.SAD.mratio.HDI.dep_inner = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 1-0.026);
% FDT.fit.SAD.mratio.HDI.dep_outer = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
% FDT.fit.SAD.mratio.HDI.dep = [FDT.fit.SAD.mratio.HDI.dep_outer(1) FDT.fit.SAD.mratio.HDI.dep_inner(2)];
% FDT.fit.SAD.mratio.HDI.dep_quantile_unc = quantile(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAD.mratio.HDIinterval_unc);
% FDT.fit.SAD.mratio.HDI.dep_inner_unc = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9);
% FDT.fit.SAD.mratio.HDI.dep_outer_unc = calc_HDI(FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
% FDT.fit.SAD.mratio.HDI.dep_unc = [FDT.fit.SAD.mratio.HDI.dep_outer_unc(1) FDT.fit.SAD.mratio.HDI.dep_inner_unc(2)];
% FDT.fit.SAD.mratio.HDI.anxVdep_quantile = quantile((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAD.mratio.HDIinterval);
% FDT.fit.SAD.mratio.HDI.anxVdep_inner = calc_HDI((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 1-0.026);
% FDT.fit.SAD.mratio.HDI.anxVdep_outer = calc_HDI((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9999);
% FDT.fit.SAD.mratio.HDI.anxVdep = [FDT.fit.SAD.mratio.HDI.anxVdep_outer(1) FDT.fit.SAD.mratio.HDI.anxVdep_inner(2)];
% FDT.fit.SAD.mratio.HDI.anxVdep_quantile_unc = quantile((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAD.mratio.HDIinterval_unc);
% FDT.fit.SAD.mratio.HDI.anxVdep_inner_unc = calc_HDI((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9);
% FDT.fit.SAD.mratio.HDI.anxVdep_outer_unc = calc_HDI((FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9999);
% FDT.fit.SAD.mratio.HDI.anxVdep_unc = [FDT.fit.SAD.mratio.HDI.anxVdep_outer_unc(1) FDT.fit.SAD.mratio.HDI.anxVdep_inner_unc(2)];
% 
% % Fit the linear model for sensitivity (filter number)
% filters_SAD = FDT.data.summary.avgFilter;
% filters_SAD(idxNaN_row_SAD) = [];
% FDT.fit.SAD.sensitivity.modelFit = fitlm(FDT.GLMs.SAD, filters_SAD);
% FDT.fit.SAD.sensitivity.coefficients = FDT.fit.SAD.sensitivity.modelFit.Coefficients.Estimate(:);
% FDT.fit.SAD.sensitivity.pvalues = FDT.fit.SAD.sensitivity.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for sensibility (confidence)
% confidence_SAD = FDT.data.summary.avgConfidence;
% confidence_SAD(idxNaN_row_SAD) = [];
% FDT.fit.SAD.sensibility.modelFit = fitlm(FDT.GLMs.SAD, confidence_SAD);
% FDT.fit.SAD.sensibility.coefficients = FDT.fit.SAD.sensibility.modelFit.Coefficients.Estimate(:);
% FDT.fit.SAD.sensibility.pvalues = FDT.fit.SAD.sensibility.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for decision bias
% FDT.fit.SAD.bias.modelFit = fitlm(FDT.GLMs.SAD, FDT.fit.SAD.mratio.modelFit.c1);
% FDT.fit.SAD.bias.coefficients = FDT.fit.SAD.bias.modelFit.Coefficients.Estimate(:);
% FDT.fit.SAD.bias.pvalues = FDT.fit.SAD.bias.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for d' (sanity check)
% FDT.fit.SAD.dPrime.modelFit = fitlm(FDT.GLMs.SAD, FDT.fit.SAD.mratio.modelFit.d1);
% FDT.fit.SAD.dPrime.coefficients = FDT.fit.SAD.dPrime.modelFit.Coefficients.Estimate(:);
% FDT.fit.SAD.dPrime.pvalues = FDT.fit.SAD.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE TAD GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the depression variable
idxDep = find(strcmp(FDT.settings.variables.summarySheet, 'DEPRESSION'));

% Create the depression variable and create GLM matrix
Dep = [];
for study = 1:4
    Dep = [Dep; FDT.data.raw.extra{study}(:,idxDep)];
end
FDT.GLMs.TAD(:,1) = TA;
FDT.GLMs.TAD(:,2) = Dep;
FDT.GLMs.TAD(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.TAD(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.TAD(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.TAD(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;

% Check for missing values and remove from GLM and FDT data
[idxNaN_row_TAD col] = find(isnan(FDT.GLMs.TAD));
FDT.GLMs.TAD(idxNaN_row_TAD,:) = [];
FDT.data.trials2counts.TAD.nR_S1 = FDT.data.trials2counts.nR_S1;
FDT.data.trials2counts.TAD.nR_S2 = FDT.data.trials2counts.nR_S2;
FDT.data.trials2counts.TAD.nR_S1(idxNaN_row_TAD) = [];
FDT.data.trials2counts.TAD.nR_S2(idxNaN_row_TAD) = [];

% Zscore variables
FDT.GLMs.TAD(:,1) = zscore(FDT.GLMs.TAD(:,1));
FDT.GLMs.TAD(:,2) = zscore(FDT.GLMs.TAD(:,2));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE TAD MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit the HMeta-d model for mratio
FDT.fit.TAD.mratio.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TAD.nR_S1, FDT.data.trials2counts.TAD.nR_S2, FDT.GLMs.TAD', FDT.settings.fit.mcmc_params);
FDT.fit.TAD.mratio.HDIinterval = [0.0 0.987];
FDT.fit.TAD.mratio.HDIinterval_unc = [0.0 0.95];
FDT.fit.TAD.mratio.tails = 1;
FDT.fit.TAD.mratio.HDI.anx_quantile = quantile(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TAD.mratio.HDIinterval);
FDT.fit.TAD.mratio.HDI.anx_inner = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
FDT.fit.TAD.mratio.HDI.anx_outer = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.TAD.mratio.HDI.anx = [FDT.fit.TAD.mratio.HDI.anx_outer(1) FDT.fit.TAD.mratio.HDI.anx_inner(2)];
FDT.fit.TAD.mratio.HDI.anx_quantile_unc = quantile(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TAD.mratio.HDIinterval_unc);
FDT.fit.TAD.mratio.HDI.anx_inner_unc = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9);
FDT.fit.TAD.mratio.HDI.anx_outer_unc = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.TAD.mratio.HDI.anx_unc = [FDT.fit.TAD.mratio.HDI.anx_outer_unc(1) FDT.fit.TAD.mratio.HDI.anx_inner_unc(2)];
FDT.fit.TAD.mratio.HDI.dep_quantile = quantile(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.TAD.mratio.HDIinterval);
FDT.fit.TAD.mratio.HDI.dep_inner = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 1-0.026);
FDT.fit.TAD.mratio.HDI.dep_outer = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.TAD.mratio.HDI.dep = [FDT.fit.TAD.mratio.HDI.dep_outer(1) FDT.fit.TAD.mratio.HDI.dep_inner(2)];
FDT.fit.TAD.mratio.HDI.dep_quantile_unc = quantile(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.TAD.mratio.HDIinterval_unc);
FDT.fit.TAD.mratio.HDI.dep_inner_unc = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9);
FDT.fit.TAD.mratio.HDI.dep_outer_unc = calc_HDI(FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.TAD.mratio.HDI.dep_unc = [FDT.fit.TAD.mratio.HDI.dep_outer_unc(1) FDT.fit.TAD.mratio.HDI.dep_inner_unc(2)];
FDT.fit.TAD.mratio.HDI.anxVdep_quantile = quantile((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAD.mratio.HDIinterval);
FDT.fit.TAD.mratio.HDI.anxVdep_inner = calc_HDI((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 1-0.026);
FDT.fit.TAD.mratio.HDI.anxVdep_outer = calc_HDI((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9999);
FDT.fit.TAD.mratio.HDI.anxVdep = [FDT.fit.TAD.mratio.HDI.anxVdep_outer(1) FDT.fit.TAD.mratio.HDI.anxVdep_inner(2)];
FDT.fit.TAD.mratio.HDI.anxVdep_quantile_unc = quantile((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAD.mratio.HDIinterval_unc);
FDT.fit.TAD.mratio.HDI.anxVdep_inner_unc = calc_HDI((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9);
FDT.fit.TAD.mratio.HDI.anxVdep_outer_unc = calc_HDI((FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAD.mratio.modelFit.mcmc.samples.mu_beta2(:)), 0.9999);
FDT.fit.TAD.mratio.HDI.anxVdep_unc = [FDT.fit.TAD.mratio.HDI.anxVdep_outer_unc(1) FDT.fit.TAD.mratio.HDI.anxVdep_inner_unc(2)];

% Fit the linear model for sensitivity (filter number)
filters_TAD = FDT.data.summary.avgFilter;
filters_TAD(idxNaN_row_TAD) = [];
FDT.fit.TAD.sensitivity.modelFit = fitlm(FDT.GLMs.TAD, filters_TAD);
FDT.fit.TAD.sensitivity.coefficients = FDT.fit.TAD.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.TAD.sensitivity.pvalues = FDT.fit.TAD.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_TAD = FDT.data.summary.avgConfidence;
confidence_TAD(idxNaN_row_TAD) = [];
FDT.fit.TAD.sensibility.modelFit = fitlm(FDT.GLMs.TAD, confidence_TAD);
FDT.fit.TAD.sensibility.coefficients = FDT.fit.TAD.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.TAD.sensibility.pvalues = FDT.fit.TAD.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.TAD.bias.modelFit = fitlm(FDT.GLMs.TAD, FDT.fit.TAD.mratio.modelFit.c1);
FDT.fit.TAD.bias.coefficients = FDT.fit.TAD.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.TAD.bias.pvalues = FDT.fit.TAD.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.TAD.dPrime.modelFit = fitlm(FDT.GLMs.TAD, FDT.fit.TAD.mratio.modelFit.d1);
FDT.fit.TAD.dPrime.coefficients = FDT.fit.TAD.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.TAD.dPrime.pvalues = FDT.fit.TAD.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE SAG GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the gender variable
idxGen = find(strcmp(FDT.settings.variables.summarySheet, 'GENDER'));

% Create the gender variable and GLM with interaction specified in terms matrix
Gender = [];
for study = 1:4
    Gender = [Gender; FDT.data.raw.extra{study}(:,idxGen)];
end
FDT.GLMs.SAG(:,1) = SA;
FDT.GLMs.SAG(:,2) = Gender;
FDT.GLMs.SAG(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.SAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;
FDT.GLMs.SAG_t = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 1 1 0 0 0];

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
FDT.fit.SAG.mratio.HDIinterval = [0.007 0.993];
FDT.fit.SAG.mratio.HDIinterval_unc = [0.025 0.975];
FDT.fit.SAG.mratio.tails = 2;
FDT.fit.SAG.mratio.SAGI.HDI.anx_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.anx = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.SAG.mratio.SAGI.HDI.anx_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.anx_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.gender_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.gender = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.987);
FDT.fit.SAG.mratio.SAGI.HDI.gender_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.gender_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG.mratio.SAGI.HDI.int_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.int = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.987);
FDT.fit.SAG.mratio.SAGI.HDI.int_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.int_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);

% Fit HDI for one-tail
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval = [0.0 0.987];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval_unc = [0.0 0.95];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_inner = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_outer = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_outer(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_inner(2)];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.9);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_unc = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_outer_unc(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.anx_inner_unc(2)];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_inner = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 1-0.026);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_outer = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_outer(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_inner(2)];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.9);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_unc = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_outer_unc(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.gender_inner_unc(2)];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_quantile = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_inner = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 1-0.026);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_outer = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_outer(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_inner(2)];
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.SAGI.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.9);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples.mu_beta3(:), 0.9999);
FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_unc = [FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_outer_unc(1) FDT.fit.SAG.mratio.SAGI.HDI.oneTail.int_inner_unc(2)];

% Fit the HMeta-d model for mratio (split interaction term)
FDT.fit.SAG.mratio.SAGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.SAG.nR_S1, FDT.data.trials2counts.SAG.nR_S2, FDT.GLMs.SAGIS', FDT.settings.fit.mcmc_params);
FDT.fit.SAG.mratio.SAGIS.HDI.anx_quantile = quantile((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.anx = calc_HDI((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.987);
FDT.fit.SAG.mratio.SAGIS.HDI.anx_quantile_unc = quantile((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.anx_unc = calc_HDI((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.gender_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.gender = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.987);
FDT.fit.SAG.mratio.SAGIS.HDI.gender_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.gender_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.int_quantile = quantile((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.int = calc_HDI((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.987);
FDT.fit.SAG.mratio.SAGIS.HDI.int_quantile_unc = quantile((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.int_unc = calc_HDI((FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.anxF_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.anxF = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.987);
FDT.fit.SAG.mratio.SAGIS.HDI.anxF_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.987);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.anxM_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.95);

% Fit HDI for one-tail
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval = [0.0 0.987];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval_unc = [0.0 0.95];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_inner = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 1-0.026);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_outer = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_outer(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_inner(2)];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.9);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta1(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_unc = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_outer_unc(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.anx_inner_unc(2)];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_inner = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 1-0.026);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_outer = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_outer(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_inner(2)];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.9);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta2(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_unc = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_outer_unc(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.gender_inner_unc(2)];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_quantile = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_inner = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 1-0.026);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_outer = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_outer(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_inner(2)];
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_quantile_unc = quantile(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.HDIinterval_unc);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_inner_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.9);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_outer_unc = calc_HDI(FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.9999);
FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_unc = [FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_outer_unc(1) FDT.fit.SAG.mratio.SAGIS.HDI.oneTail.int_inner_unc(2)];

% Fit the linear model for sensitivity (filter number)
filters_SAG = FDT.data.summary.avgFilter;
filters_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG.sensitivity.modelFit = fitlm(FDT.GLMs.SAGI, filters_SAG, 'CategoricalVars', 2);
FDT.fit.SAG.sensitivity.coefficients = FDT.fit.SAG.sensitivity.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensitivity.pvalues = FDT.fit.SAG.sensitivity.modelFit.Coefficients.pValue(:);

% Fit the linear model for sensibility (confidence)
confidence_SAG = FDT.data.summary.avgConfidence;
confidence_SAG(idxNaN_row_SAG) = [];
FDT.fit.SAG.sensibility.modelFit = fitlm(FDT.GLMs.SAGI, confidence_SAG, 'CategoricalVars', 2);
FDT.fit.SAG.sensibility.coefficients = FDT.fit.SAG.sensibility.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.sensibility.pvalues = FDT.fit.SAG.sensibility.modelFit.Coefficients.pValue(:);

% Fit the linear model for decision bias
FDT.fit.SAG.bias.modelFit = fitlm(FDT.GLMs.SAGI, FDT.fit.SAG.mratio.SAGI.modelFit.c1, 'CategoricalVars', 2);
FDT.fit.SAG.bias.coefficients = FDT.fit.SAG.bias.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.bias.pvalues = FDT.fit.SAG.bias.modelFit.Coefficients.pValue(:);

% Fit the linear model for d' (sanity check)
FDT.fit.SAG.dPrime.modelFit = fitlm(FDT.GLMs.SAGI, FDT.fit.SAG.mratio.SAGI.modelFit.d1, 'CategoricalVars', 2);
FDT.fit.SAG.dPrime.coefficients = FDT.fit.SAG.dPrime.modelFit.Coefficients.Estimate(:);
FDT.fit.SAG.dPrime.pvalues = FDT.fit.SAG.dPrime.modelFit.Coefficients.pValue(:);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE TAG GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the gender GLM with interaction specified in terms matrix
FDT.GLMs.TAG(:,1) = TA;
FDT.GLMs.TAG(:,2) = Gender;
FDT.GLMs.TAG(1:length(FDT.data.IDs{1}),3:5) = -1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2}),3) = 1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3}),4) = 1;
FDT.GLMs.TAG(length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+1:length(FDT.data.IDs{1})+length(FDT.data.IDs{2})+length(FDT.data.IDs{3})+length(FDT.data.IDs{4}),5) = 1;
FDT.GLMs.TAG_t = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 1 1 0 0 0];

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

% % Fit the HMeta-d model for mratio (specified interaction term)
% FDT.fit.TAG.mratio.TAGI.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TAG.nR_S1, FDT.data.trials2counts.TAG.nR_S2, FDT.GLMs.TAGI', FDT.settings.fit.mcmc_params);
% FDT.fit.TAG.mratio.HDIinterval = [0.007 0.993];
% FDT.fit.TAG.mratio.HDIinterval_unc = [0.025 0.975];
% FDT.fit.TAG.mratio.tails = 2;
% FDT.fit.TAG.mratio.TAGI.HDI.anx_quantile = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGI.HDI.anx = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:), 0.987);
% FDT.fit.TAG.mratio.TAGI.HDI.anx_quantile_unc = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGI.HDI.anx_unc = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta1(:), 0.95);
% FDT.fit.TAG.mratio.TAGI.HDI.gender_quantile = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGI.HDI.gender = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:), 0.987);
% FDT.fit.TAG.mratio.TAGI.HDI.gender_quantile_unc = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGI.HDI.gender_unc = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta2(:), 0.95);
% FDT.fit.TAG.mratio.TAGI.HDI.int_quantile = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGI.HDI.int = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:), 0.987);
% FDT.fit.TAG.mratio.TAGI.HDI.int_quantile_unc = quantile(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGI.HDI.int_unc = calc_HDI(FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples.mu_beta3(:), 0.95);
% 
% % Fit the HMeta-d model for mratio (split interaction term)
% FDT.fit.TAG.mratio.TAGIS.modelFit = fit_meta_d_mcmc_regression(FDT.data.trials2counts.TAG.nR_S1, FDT.data.trials2counts.TAG.nR_S2, FDT.GLMs.TAGIS', FDT.settings.fit.mcmc_params);
% FDT.fit.TAG.mratio.TAGIS.HDI.anx_quantile = quantile((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGIS.HDI.anx = calc_HDI((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.987);
% FDT.fit.TAG.mratio.TAGIS.HDI.anx_quantile_unc = quantile((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGIS.HDI.anx_unc = calc_HDI((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) + FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
% FDT.fit.TAG.mratio.TAGIS.HDI.gender_quantile = quantile(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGIS.HDI.gender = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.987);
% FDT.fit.TAG.mratio.TAGIS.HDI.gender_quantile_unc = quantile(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta3(:), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGIS.HDI.gender_unc = calc_HDI(FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta3(:), 0.95);
% FDT.fit.TAG.mratio.TAGIS.HDI.int_quantile = quantile((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAG.mratio.HDIinterval);
% FDT.fit.TAG.mratio.TAGIS.HDI.int = calc_HDI((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.987);
% FDT.fit.TAG.mratio.TAGIS.HDI.int_quantile_unc = quantile((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), FDT.fit.TAG.mratio.HDIinterval_unc);
% FDT.fit.TAG.mratio.TAGIS.HDI.int_unc = calc_HDI((FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta1(:) - FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples.mu_beta2(:)), 0.95);
% 
% % Fit the linear model for sensitivity (filter number)
% filters_TAG = FDT.data.summary.avgFilter;
% filters_TAG(idxNaN_row_TAG) = [];
% FDT.fit.TAG.sensitivity.modelFit = fitlm(FDT.GLMs.TAGI, filters_TAG, 'CategoricalVars', 2);
% FDT.fit.TAG.sensitivity.coefficients = FDT.fit.TAG.sensitivity.modelFit.Coefficients.Estimate(:);
% FDT.fit.TAG.sensitivity.pvalues = FDT.fit.TAG.sensitivity.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for sensibility (confidence)
% confidence_TAG = FDT.data.summary.avgFilter;
% confidence_TAG(idxNaN_row_TAG) = [];
% FDT.fit.TAG.sensibility.modelFit = fitlm(FDT.GLMs.TAGI, confidence_TAG, 'CategoricalVars', 2);
% FDT.fit.TAG.sensibility.coefficients = FDT.fit.TAG.sensibility.modelFit.Coefficients.Estimate(:);
% FDT.fit.TAG.sensibility.pvalues = FDT.fit.TAG.sensibility.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for decision bias
% FDT.fit.TAG.bias.modelFit = fitlm(FDT.GLMs.TAGI, FDT.fit.TAG.mratio.TAGIS.modelFit.c1, 'CategoricalVars', 2);
% FDT.fit.TAG.bias.coefficients = FDT.fit.TAG.bias.modelFit.Coefficients.Estimate(:);
% FDT.fit.TAG.bias.pvalues = FDT.fit.TAG.bias.modelFit.Coefficients.pValue(:);
% 
% % Fit the linear model for d' (sanity check)
% FDT.fit.TAG.dPrime.modelFit = fitlm(FDT.GLMs.TAGI, FDT.fit.TAG.mratio.TAGIS.modelFit.d1, 'CategoricalVars', 2);
% FDT.fit.TAG.dPrime.coefficients = FDT.fit.TAG.dPrime.modelFit.Coefficients.Estimate(:);
% FDT.fit.TAG.dPrime.pvalues = FDT.fit.TAG.dPrime.modelFit.Coefficients.pValue(:);


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
% FDT.fit.Dep.mratio.modelFit.mcmc.samples = [];
% FDT.fit.SAD.mratio.modelFit.mcmc.samples = [];
% FDT.fit.TAD.mratio.modelFit.mcmc.samples = [];
% FDT.fit.SAG.mratio.SAGI.modelFit.mcmc.samples = [];
% FDT.fit.SAG.mratio.SAGIS.modelFit.mcmc.samples = [];
% FDT.fit.TAG.mratio.TAGI.modelFit.mcmc.samples = [];
% FDT.fit.TAG.mratio.TAGIS.modelFit.mcmc.samples = [];

% Save data matrix
save(FDT.settings.names.output, 'FDT', '-v7.3');


end

