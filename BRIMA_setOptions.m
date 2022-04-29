%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SET OPTIONS FOR BRIMA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Olivia Harrison
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
addpath(genpath('/Users/faullol/Documents/just_breathe/otago/OneDrive - University of Otago/code/BRIMA'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY DATA PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify data path
options.paths.FDTdata = fullfile(options.paths.root, 'BRIMA-FDTdata');

% Specify output path for analyses
options.paths.analysis = fullfile(options.paths.root, 'BRIMA-results');

% Specify output path for figures
options.paths.figures = fullfile(options.paths.root, 'BRIMA-figures');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY DATA NAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify combined EXCEL data
options.names.dataFile = fullfile(options.paths.root, 'BRIMA-scores.xlsx');

% Specify data sheets for each study
options.names.sheetNames = {'BATH', 'BIRMINGHAM', 'ZURICH', 'OXFORD'};

% Specify data columns
options.variables.summarySheet = {'AGE', 'GENDER', 'ANXIETY-S', 'ANXIETY-T', 'DEPRESSION'};

% Specify data analysis output
options.names.output = fullfile(options.paths.analysis, 'BRIMA-output.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRA SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number of confidence bins
options.confidenceBinEdges = linspace(0,100,11);
options.confidenceBins = 10;

% Specify H-Medtad parameters
options.fit.mcmc_params = fit_meta_d_params;
options.fit.mcmc_params.estimate_dprime = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end