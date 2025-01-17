function LEiDA_AnalysisK_JL(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse the PL states chosen according to the analysis from
% LEiDA_Start. The user should select a value for K and this function will
% return summary plots of the set of PL states.
%
% This function contains two sections: 
%       (A) User defines the parameters and selects the value of K.
%       (B) Generate and save figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis plots for selected K:
%    - Plot PL states as links in cortex
%    - Plot PL states as nodes in cortex
%    - Plot PL states in 3D glass brain
%    - Plot PL states in matrix format
%    - Plot PL states in vector format with number of each parcel
%    - Plot boxplot of fractional occupancy values
%    - Plot boxplot of dwell time values
%    - Plot of summary information for the repertoire of K PL states
%    - Plot overlap of PL states with RSNs defined by Yeo et al., (2011)
%
% Tutorial: README.md
% Version:  V1.0, June 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

%% A: USER INPUT PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = param.cluster.SelectK;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = param.LEiDA_directory;
% Name of the run to be used to create the folder to save the data:
res = fullfile(param.res.atlas, param.res.preproc, param.res.run_name);

% Parcellation used to run LEiDA_Start script:
Parcellation = param.parcellation;

% Add the LEiDA_directory to the matlab path
addpath(genpath(LEiDA_directory))

%% B: ANALYSE K CENTROIDS SELECTED ACCORDING TO OUTPUT FROM LEiDA_Start

% Close all open figures
close all;

% Go to the directory containing the LEiDA functions
%cd(LEiDA_directory) %- no need? has been added to path JL

% Directory with the results from LEiDA
%leida_res = fullfile(LEiDA_directory, ['res_' run_name]); % changed by JL to make it more general and not require specification of '/'
[pathstr,~,~] = fileparts(LEiDA_directory); % added by JL
mainPath = pathstr(1:find(pathstr == filesep, 1, 'last')); % added by JL
leida_res = fullfile(mainPath, 'results',res); % added by JL

% Create a directory to store results for defined value of K
K_dir = fullfile(leida_res, ['K' num2str(SelectK)]); % changed by JL to make it more general and not require specification of '/'

if ~exist(K_dir, 'dir')
    mkdir(K_dir);
end


disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR K = ' num2str(SelectK) ' PL STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Plot PL states as links in cortex
Plot_K_links_in_cortex(leida_res,K_dir,SelectK,Parcellation);

% Plot PL states as nodes in cortex
Plot_K_nodes_in_cortex(leida_res,K_dir,SelectK,Parcellation);

% Plot PL states in 3D glass brain
Plot_K_3Dbrain(leida_res,K_dir,SelectK,Parcellation);

% Plot PL states in matrix format
Plot_K_matrix(leida_res,K_dir,SelectK);

% Plot PL states in vector format with number of each parcel
Plot_K_vector_numbered(leida_res,K_dir,SelectK);

% Plot boxplot of fractional occupancy values
Plot_K_boxplot_FO(leida_res,K_dir,SelectK,Parcellation);

% Plot boxplot of dwell time values
Plot_K_boxplot_DT(leida_res,K_dir,SelectK,Parcellation);

% Plot of summary information for the repertoire of K PL states
Plot_K_repertoire(leida_res,K_dir,SelectK,Parcellation);

% Plot overlap between PL states and Yeo RSNs
Plot_K_overlap_yeo_nets(leida_res,K_dir,SelectK,Parcellation)
