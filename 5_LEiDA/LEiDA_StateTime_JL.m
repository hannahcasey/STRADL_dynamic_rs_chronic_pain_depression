function LEiDA_StateTime_JL(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to plot the state time courses for all participants and plot the
% state time courses for a specific subject. The user should select a value
% for K and a subject to analyse.
%
% This function contains two sections:
%       (A) User defines the parameters and selects the value of K.
%       (B) Generate and save figures of the state time curses for the
%           selected K and subject.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis plots for selected K and subject:
%    - Plot state time courses for all participants
%    - Plot cluster blocks on the fMRI time series for specific subject
%    - Plot cluster stairs for specific subject
%
% Tutorial: README.md
% Version:  V1.0, June 2022
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

%% A: USER INPUT PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = param.cluster.SelectK;
% Define the subject to analyse into more detail (file name of subject or unique ID/number):
Subject = param.Subject;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = param.LEiDA_directory;
% Name of the run to be used to create the folder to save the data:
res = fullfile(param.res.atlas, param.res.preproc, param.res.run_name);


% Add the LEiDA_directory to the matlab path
addpath(genpath(LEiDA_directory))

%% B: PLOT STATE TIME COURSES FOR K CENTROIDS AND/OR A SPECIFIC SUBJECT
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
K_dir = fullfile(leida_res, ['K' num2str(SelectK)]);

if ~exist(K_dir, 'dir')
    mkdir(K_dir);
end

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATE TIME COURSES FOR K = ' num2str(SelectK) ' CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Plot state time courses for all subjects by condition
Plot_K_state_time(leida_res,K_dir,SelectK);

% Plot state time course for specific subject and fMRI signal
Plot_subj_cluster_blocks(leida_res,K_dir,SelectK,Subject);

% Plot state time course for specific subject as stairs plot
Plot_subj_stairs(leida_res,K_dir,SelectK,Subject);

