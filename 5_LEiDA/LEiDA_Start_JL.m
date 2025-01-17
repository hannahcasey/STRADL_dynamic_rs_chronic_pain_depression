function LEiDA_Start_JL(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function runs LEiDA on any given dataset.
%
% Since we do not know a priori the optimal number of PL states that
% differentiate between conditions, this function analyses the solutions
% obtained for different numbers of PL states (K).
%
% After analysing the output figures saved in the folder LEiDA_Results, the
% user can then choose the optimal number of PL states for subsequent
% detailed analyses using the functions LEiDA_AnalysisK.m and
% LEiDA_AnalysisCentroid.m.
%
% This function contains two sections:
%       (A) User defines the parameters and properties of the study.
%       (B) Run LEiDA and statistics.
%       (C) Generate and save figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis:
%    - Compute the leading eigenvectors for all participants
%    - Cluster the leading eigenvectors of all participants
%    - Compute statistics to compare across conditions
% C: Figures
%    - Analysis of Fractional Occupancy values
%    - Analysis of Dwell Time values
%    - Pyramid of PL states
%
% Tutorial: README.md
% Version:  V1.0, June 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

%% A: STUDY PARAMETERS
% section A changed by JL to receive input of function instead

% Directory of the LEiDA toolbox folder:
LEiDA_directory = param.LEiDA_directory;
% Directory of the folder with the parcellated neuroimaging data:
Data_directory = param.Data_directory;
% Name of the run to be used to create the folder to save the data:
res = fullfile(param.res.atlas, param.res.preproc, param.res.run_name);

% Tag of conditions given in the parcellated image files:
Conditions_tag = param.Conditions_tag; %'placebo_pre','placebo_post','menthol_pre','menthol_post'}; 
                                                % choose 2 of the 4, select paired or unpaired accordingly
% Parcellation applied to the imaging data (see tutorial):
Parcellation = param.parcellation;

% Number of brain areas to consider for analysis:
N_areas = param.parcellation.n_areas; % added by LK 15/02/2023

% Repetition time (TR) of the fMRI data (if unknown set to 1):
TR = param.fMRI.TR;
% Maximum number of TRs for all fMRI sessions:
Tmax = param.fMRI.Tmax;
% Apply temporal filtering to data (0: no; 1: yes)
apply_filter = param.fMRI.apply_filter;
% Lowpass frequency of filter (default 0.1):
flp = param.fMRI.flp;
% Highpass frequency of filter (default 0.01):
fhi = param.fMRI.fhi;

% For the statistics:
% Choose 0 (unpaired) if subjects in different conditions are not the
% same; or 1 (paired) if subjects are the same across conditions.
Paired_tests = param.stats.Paired_tests;
% Number of permutations. For the first analysis to be relatively quick,
% run around 500 permutations, but then increase to 10000 to increase the
% reliability of the final statistical results (p-values) for publication.
n_permutations = param.stats.n_permutations;
% Number of bootstrap samples within each permutation. For the first
% analysis to be relatively quick, choose around 10, but then increase to
% 500 for more reliable final results.
n_bootstraps = param.stats.n_bootstraps;

% For the figure of the pyramid of PL states:
% Direction to plot the FC states/brain ('SideView' or 'TopView'):
CortexDirection = param.plot.CortexDirection;


% Add the LEiDA_directory to the matlab path
addpath(genpath(LEiDA_directory))

%% B: RUN LEADING EIGENVECTOR DYNAMICS ANALYSIS

% Go to the directory containing the LEiDA functions
%cd(LEiDA_directory) %- no need? has been added to path JL

% Create a directory to store the results from the current LEiDA run
%leida_res = fullfile(LEiDA_directory, ['res_' run_name]); % changed by JL to make it more general and not require specification of '/'
[pathstr,~,~] = fileparts(LEiDA_directory); % added by JL
mainPath = pathstr(1:find(pathstr == filesep, 1, 'last')); % added by JL
leida_res = fullfile(mainPath, 'results',res); % added by JL

if ~exist(leida_res, 'dir') % changed by JL to make it more general and not require specification of '/'
    mkdir(leida_res); % changed by JL to make it more general and not require specification of '/'
end



% Compute the leading eigenvectors of the data
if exist(fullfile(leida_res, 'LEiDA_EigenVectors.mat'),'file') % added by JL so to run only if t hasn't yet or if you want to run again
    run_again = input('LEiDA_data has already been run. Do you want to run again? (y/n)','s');
    if strcmp(run_again,'y')
        LEiDA_data(Data_directory,leida_res,N_areas,Tmax,apply_filter,flp,fhi,TR);
    end
else
    disp('running LEIDA data')
    tic
    LEiDA_data(Data_directory,leida_res,N_areas,Tmax,apply_filter,flp,fhi,TR);
    toc
    disp('finished LEIDA data')
end

% Cluster the leading eigenvectors of all subjects
if exist(fullfile(leida_res, 'LEiDA_Clusters.mat'),'file') % added by JL so to run only if t hasn't yet or if you want to run again
    run_again = input('LEiDA_cluster has already been run. Do you want to run again? (y/n)','s');
    if strcmp(run_again,'y')
        LEiDA_cluster(leida_res);
    end
else
    disp('running LEIDA data')
    tic
    LEiDA_cluster(leida_res);
    toc
end

% Compute the fractional occupancy and perform hypothesis tests
if exist(fullfile(leida_res, 'LEiDA_Stats_FracOccup.mat'),'file') % added by JL so to run only if t hasn't yet or if you want to run again
    run_again = input('LEiDA_stats_FracOccup has already been run. Do you want to run again? (y/n)','s');
    if strcmp(run_again,'y')
        LEiDA_stats_FracOccup(leida_res,Conditions_tag,Paired_tests,n_permutations,n_bootstraps);
    end
else
    LEiDA_stats_FracOccup(leida_res,Conditions_tag,Paired_tests,n_permutations,n_bootstraps);
end

% Compute the dwell time and perform hypothesis tests
if exist(fullfile(leida_res, 'LEiDA_Stats_DwellTime.mat'),'file') % added by JL so to run only if t hasn't yet or if you want to run again
    run_again = input('LEiDA_stats_DwellTime has already been run. Do you want to run again? (y/n)','s');
    if strcmp(run_again,'y')
        LEiDA_stats_DwellTime(leida_res,Conditions_tag,Paired_tests,TR,n_permutations,n_bootstraps);
    end
else
    LEiDA_stats_DwellTime(leida_res,Conditions_tag,Paired_tests,TR,n_permutations,n_bootstraps);
end

%% C: MAKE FIGURES

% Generate and save the p-value and barplot plots for fractional occupancy
Plot_FracOccup(leida_res)

% Generate and save the p-value and barplot plots for dwell time
Plot_DwellTime(leida_res)

% Plot the centroids obtained using LEiDA and their overlap with Yeo nets
Plot_Centroid_Pyramid(leida_res,Conditions_tag,Parcellation,CortexDirection)