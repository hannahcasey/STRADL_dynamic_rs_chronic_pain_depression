function LEiDA_TransitionsK_JL(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the transition probability matrix for each subject
% and compare the transition probabilities between conditions for a value
% of K chosen according to the analysis from LEiDA_Start.
%
% This function contains two sections:
%       (A) User defines the parameters and properties of the study.
%       (B) Compute the transition probability matrix and intergroup 
%           differences in transition probabilities.
%       (C) Generate and save figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis for selected K:
%    - Compute Transition Probability Matrix (TPM) for each participant
%    - Compare state-to-state transition probabilities between conditions
% C: Figures
%    - Plot the mean TPM for each condition
%    - Plot the intergroup differences between transition probabilties
%
% Tutorial: README.md
% Version:  V1.0, June 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

%% A: STUDY PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = param.cluster.SelectK;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = param.LEiDA_directory;
% Name of the run to be used to create the folder to save the data:
res = fullfile(param.res.atlas, param.res.preproc, param.res.run_name);

% For the statistics:
% Number of permutations. For the first analysis to be relatively quick,
% run around 500 permutations, but then increase to 10000 to increase the
% reliability of the final statistical results (p-values) for publication.
n_permutations = param.stats.n_permutations;

% Number of bootstrap samples within each permutation. For the first
% analysis to be relatively quick, choose around 10, but then increase to
% 500 for more reliable final results.
n_bootstraps = param.stats.n_bootstraps;

% Add the LEiDA_directory to the matlab path
addpath(genpath(LEiDA_directory))

%% B: COMPUTE TPM AND ANALYSE INTERGROUP DIFFERENCES IN TRANSITION PROBABILITIES FOR SELECTED K

% Close all open figures
close all;

% Go to the directory containing the LEiDA functions
%cd(LEiDA_directory) %- no need? has been added to path JL

% Directory with the results from LEiDA
%leida_res = fullfile(LEiDA_directory, ['res_' run_name]); % changed by JL to make it more general and not require specification of '/'
[pathstr,~,~] = fileparts(LEiDA_directory); % added by JL
mainPath = pathstr(1:find(pathstr == filesep, 1, 'last')); % added by JL
leida_res = fullfile(mainPath, 'results', res); % added by JL


% Create a directory to store results for defined value of K
K_dir = fullfile(leida_res, ['K' num2str(SelectK)]);  % changed by JL to make it more general

if ~exist(K_dir, 'dir')
    mkdir(K_dir);
end

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITIONS FOR K = ' num2str(SelectK) ' CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Compute the transition probability matrix and perform hypothesis tests
LEiDA_stats_TransitionMatrix(leida_res,K_dir,SelectK,n_permutations,n_bootstraps);

%% C: MAKE FIGURES

% Plot the mean transition probability matrix
Plot_K_tpm(K_dir,SelectK);

% Plot summary of differences in transition probabilities between conditions
Plot_K_diffs_transitions(K_dir,SelectK);