%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SET OPTIONS FOR BRIMA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Olivia Faull
% Created: 12/04/2022
% -------------------------------------------------------------------------
% TO RUN:   options     = PBIHB_setOptions
% OUTPUTS:  options     = Structure with all specified options saved for
%                         use in further analysis
% -------------------------------------------------------------------------
% DESCRIPTION:
% This script sets the required options, paths and file names for the FDT 
% analysis for the BRIMA study (combined study sites). Any required changes
% should be made here.
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function options = BRIMA_setOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY DATA PATH -- TO ADJUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify path root for data
options.paths.root = '/Volumes/sci-psychology/Administration/Olivia_Harrison/study-BRIMA';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD REQUIRED CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add BRIMA code folder
addpath(genpath('../HMeta-d/'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY DATA PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify data path
options.paths.FDTdata = fullfile(options.paths.root, 'FDTdata');

% Specify output path for analyses
options.paths.analysis = fullfile(options.paths.root, 'BRIMA-results');

% Specify output path for figures
options.paths.figures = fullfile(options.paths.root, 'BRIMA-figures');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY OPTIONS FOR PCA vs EFA SCORES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.pca.asthmaPCA = 0; % Change to 1 if using PCA scores for asthma
options.pca.asthmaEFA = 1; % Change to 1 if using EFA scores for asthma


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY NAMES AND SAVE NAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify data files
options.names.fullDataFile = fullfile(options.paths.data, 'Asthma_imputed_2021.xlsx');
options.names.sheetAsth = 'Asthma_imputed';
options.names.sheetCont = 'HV_imputedFiles';
options.names.sheetAsthRaw = 'Asthma_rawFiles';
options.names.sheetContRaw = 'HV_rawFiles';
% options.names.sheetAsthPPID = 'AS_IDS';
% options.names.sheetContPPID = 'HV_IDS';
% options.names.sheetAsthQ = 'AS_imputedFiles';
% options.names.sheetAsthQRaw = 'Asthma_rawFiles';
% options.names.sheetContQ = 'HV_imputedFiles';
% options.names.sheetContQRaw = 'HV_rawFiles';
% options.names.sheetAsthExtra = 'AS_Extra';
% options.names.sheetAsthExtraRaw = 'AS_rawExtra';
% options.names.sheetContExtra = 'HV_Extra';
% options.names.sheetContExtraRaw = 'HV_rawExtra';
options.names.sheetAsthEFA = 'AS_EFA';
options.names.FDTnameEnd = '_Interoception.xlsx';

% Specify questionnaire variable names
options.names.questionnaires.names.short = {'stai', 'trai', 'cesd', 'asi',...
    'hai', 'maiaAtt', 'maiaBL', 'maiaEA', 'maiaND', 'maiaNot', 'maiaNW',...
    'maiaSR', 'maiaTr', 'nijm', 'mas', 'd12', 'fss', 'bmqCon', 'bmqHarm', 'bmqNec',...
    'bmqOU', 'caaEx', 'caaGen', 'acqCon', 'acqlqEm', 'acqlqEn', 'acqlqLim',...
    'acqlqSymp'};
options.names.questionnaires.names.long = {'quest.stai',...
    'quest.trai', 'quest.cesd', 'quest.asi',...
    'quest.hai', 'quest.maiaAttentionRegulation',...
    'quest.maiaBodyListening', 'quest.maiaEmotionalAwareness',...
    'quest.maiaNotDistracting', 'quest.maiaNoticing',...
    'quest.maiaNotWorrying', 'quest.maiaSelfRegulation',...
    'quest.maiaTrusting', 'quest.nijmegan',...
    'quest.masMedication', 'quest.d12',...
    'quest.fssFatigue', 'quest.bmqConcerns',...
    'quest.bmqHarm', 'quest.bmqNecessity',...
    'quest.bmqOveruse', 'quest.caaExacerbation',...
    'quest.caaGeneral', 'quest.acqControl',...
    'quest.acqlqEmotional', 'quest.acqlqEnvironmental',...
    'quest.acqlqLimitations', 'quest.acqlqSymptoms'};
options.names.questionnaires.names.HVlong = {'questionnaires.staiTotal',...
    'quest.trai', 'quest.cesd', 'quest.asi',...
    'quest.hai', 'quest.maiaAttentionRegulation',...
    'quest.maiaBodyListening', 'quest.maiaEmotionalAwareness',...
    'quest.maiaNotDistracting', 'quest.maiaNoticing',...
    'quest.maiaNotWorrying', 'quest.maiaSelfRegulation',...
    'quest.maiaTrusting', 'quest.nijmegan',...
    'quest.masMedication', 'quest.d12', 'quest.fssFatigue'};
options.names.questionnaires.names.stai = 'quest.stai';
options.names.questionnaires.names.trai = 'quest.trai';
options.names.questionnaires.names.cesd = 'quest.cesd';
options.names.questionnaires.names.asi = 'quest.asi';
options.names.questionnaires.names.hai = 'quest.hai';
options.names.questionnaires.names.maiaAtt = 'quest.maiaAttentionRegulation';
options.names.questionnaires.names.maiaBL = 'quest.maiaBodyListening';
options.names.questionnaires.names.maiaEA = 'quest.maiaEmotionalAwareness';
options.names.questionnaires.names.maiaND = 'quest.maiaNotDistracting';
options.names.questionnaires.names.maiaNot = 'quest.maiaNoticing';
options.names.questionnaires.names.maiaNW = 'quest.maiaNotWorrying';
options.names.questionnaires.names.maiaSR = 'quest.maiaSelfRegulation';
options.names.questionnaires.names.maiaTr = 'quest.maiaTrusting';
options.names.questionnaires.names.nijm = 'quest.nijmegan';
options.names.questionnaires.names.mas = 'quest.masMedication';
options.names.questionnaires.names.d12 = 'quest.d12';
options.names.questionnaires.names.bmqCon = 'quest.bmqConcerns';
options.names.questionnaires.names.bmqHarm = 'quest.bmqHarm';
options.names.questionnaires.names.bmqNec = 'quest.bmqNecessity';
options.names.questionnaires.names.bmqOU = 'quest.bmqOveruse';
options.names.questionnaires.names.caaEx = 'quest.caaExacerbation';
options.names.questionnaires.names.caaGen = 'quest.caaGeneral';
options.names.questionnaires.names.acqCon = 'quest.acqControl';
options.names.questionnaires.names.acqlqEm = 'quest.acqlqEmotional';
options.names.questionnaires.names.acqlqEn = 'quest.acqlqEnvironmental';
options.names.questionnaires.names.acqlqLim = 'quest.acqlqLimitations';
options.names.questionnaires.names.acqlqSymp = 'quest.acqlqSymptoms';
options.names.questionnaires.names.fss = 'quest.fssFatigue';

% Specify variables to be reverse-scored
options.names.questionnaires.reversed = {options.names.questionnaires.names.mas, options.names.questionnaires.names.acqCon, options.names.questionnaires.names.acqlqEm, options.names.questionnaires.names.acqlqEn, options.names.questionnaires.names.acqlqLim, options.names.questionnaires.names.acqlqSymp};

% Specify physiological variable names
options.names.phys.long = {'screening.physicalFEV1', 'screening.physicalFVC', ...
    'bronchodilation.afterFEV1', 'bronchodilation.beforeFEV1', 'bronchodilation.afterFVC', ...
    'bronchodilation.beforeFVC', 'screening.physicalFEV1PercentPredicted', ...
    'screening.physicalPEFR', 'exhaled.exhaledNoFeNo', 'bloods.bloodsEosinCount', ...
    'screening.physicalHeight', 'screening.physicalWeight', 'general.generalGender', ...
    'general.generalAge'};
options.names.phys.FEV1 = 'screening.physicalFEV1';
options.names.phys.FVC = 'screening.physicalFVC';
options.names.phys.broncho.afterFEV1 = 'bronchodilation.afterFEV1';
options.names.phys.broncho.beforeFEV1 = 'bronchodilation.beforeFEV1';
options.names.phys.broncho.afterFVC = 'bronchodilation.afterFVC';
options.names.phys.broncho.beforeFVC = 'bronchodilation.beforeFVC';
options.names.phys.FEV1Pred = 'screening.physicalFEV1PercentPredicted';
options.names.phys.peakFlow = 'screening.physicalPEFR';
options.names.phys.FeNO = 'exhaled.exhaledNoFeNo';
options.names.phys.bloodEos = 'bloods.bloodsEosinCount';
options.names.phys.height = 'screening.physicalHeight';
options.names.phys.weight = 'screening.physicalWeight';
options.names.phys.gender = 'general.generalGender';
options.names.phys.age = 'general.generalAge';

% Specify attention task variable names
options.names.attention.long = {'alerting', 'orienting', 'executive'};
options.names.attention.alerting = 'alerting';
options.names.attention.orienting = 'orienting';
options.names.attention.executive = 'executive';
options.names.VPT.long = {'asthmaBias'};
options.names.attention.asthmaBias = 'asthmaBias';
options.names.attention.asthmaVigilance = 'asthmaVigilance';
options.names.attention.asthmaEngage = 'asthmaEngage';

% Specify save name for data analysis
options.saveNames.analysis.asthmaCluster = fullfile(options.paths.analysis, 'asthma_cluster_analysis');
options.saveNames.analysis.controlsCluster = fullfile(options.paths.analysis, 'controls_cluster_analysis');
options.saveNames.analysis.allPCAs = fullfile(options.paths.analysis, 'total_PCA');
options.saveNames.analysis.asthmaStrat = fullfile(options.paths.analysis, 'asthma_stratification_analysis');
options.saveNames.analysis.asthmaFDT = fullfile(options.paths.analysis, 'asthma_FDT_analysis');
options.saveNames.analysis.regression = fullfile(options.paths.analysis, 'total_regression_analysis');
options.saveNames.analysis.logRegression = fullfile(options.paths.analysis, 'asthma_logRegression_analysis');
options.saveNames.analysis.groupDiffs = fullfile(options.paths.analysis, 'asthma_groupDiffs_analysis');
options.saveNames.analysis.asthmaClusterStratFull = fullfile(options.paths.analysis, 'asthma_cluster_analysis_full');

% Specify save names for raw figures
options.saveNames.figures.asthmaClustergram_measures = fullfile(options.paths.figuresRaw, 'asthma_clustergram_measures');
options.saveNames.figures.controlsClustergram_measures = fullfile(options.paths.figuresRaw, 'controls_clustergram_measures');
options.saveNames.figures.asthmaClustergram_subjects = fullfile(options.paths.figuresRaw, 'asthma_clustergram_subjects');
options.saveNames.figures.asthmaClustergram_subjectsElbow = fullfile(options.paths.figuresRaw, 'asthma_clustergram_subjectsElbow');
options.saveNames.figures.asthmaClustergram_measuresTotal = fullfile(options.paths.figuresRaw, 'asthma_clustergram_measuresTotal');
options.saveNames.figures.asthmaClustergram_subjectsTotal = fullfile(options.paths.figuresRaw, 'asthma_clustergram_subjectsTotal');

% Specify names for final figures
options.saveNames.figures.metacog.asthmaFactors = fullfile(options.paths.figures, 'asthma_FDT_figure_metacog');
options.saveNames.figures.compare.totalFactors = fullfile(options.paths.figures, 'asthmaVcontrols_figure_regressions');
options.saveNames.figures.compare.asthmaControls.PCA = fullfile(options.paths.figures, 'asthmaVcontrols_figure_PCA');
options.saveNames.figures.compare.asthmaControls.phys = fullfile(options.paths.figures, 'asthmaVcontrols_figure_phys');
options.saveNames.figures.compare.asthmaControls.extra = fullfile(options.paths.figures, 'asthmaVcontrols_figure_extra');
options.saveNames.figures.compare.asthmaGroups.PCA = fullfile(options.paths.figures, 'asthmaGroups_figure_PCA');
options.saveNames.figures.compare.asthmaGroups.phys = fullfile(options.paths.figures, 'asthmaGroups_figure_phys');
options.saveNames.figures.compare.asthmaGroups.extra = fullfile(options.paths.figures, 'asthmaGroups_figure_extra');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY FIGURE DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.figDim.width1 = 500;
options.figDim.width2 = 670;
options.figDim.width3 = 1010;
options.figDim.width4 = 1350;
options.figDim.width5 = 1500;
options.figDim.height1 = 370;
options.figDim.height2 = 630;
options.figDim.height3 = 840;
options.figDim.height4 = 1050;
options.figDim.height5 = 1350;


end