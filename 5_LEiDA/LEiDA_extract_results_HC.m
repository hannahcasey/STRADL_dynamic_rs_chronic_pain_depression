function LEiDA_extract_results_HC(param)

%% This function extracts and saves all results needed form LEiDA run

%% Get file paths
% Directory of the LEiDA toolbox folder:
LEiDA_directory = param.LEiDA_directory;
% Directory of the folder with the parcellated neuroimaging data:
Data_directory = param.Data_directory;
% Name of the run to be used to create the folder to save the data:
res = fullfile(param.res.atlas, param.res.preproc, param.res.run_name);

[pathstr,~,~] = fileparts(LEiDA_directory); % added by JL
mainPath = pathstr(1:find(pathstr == filesep, 1, 'last')); % added by JL
leida_res = fullfile(mainPath, 'results', res); % added by JL

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = param.cluster.SelectK;

% Create a directory to store results for defined value of K
K_dir = fullfile(leida_res, ['K' num2str(SelectK)]);

% Make a new results folder
extract_results_dir = fullfile(mainPath, 'results/extracted_results');

if ~exist(extract_results_dir, 'dir') 
    mkdir(extract_results_dir); 
end

%% Copy results to new results folder
% Cluster plots
OverlapYeoNets = imread(fullfile([K_dir '/K' num2str(SelectK), '_OverlapYeoNets.png']));
imwrite(OverlapYeoNets, fullfile(extract_results_dir, 'OverlapYeoNets.png'));

cluster_brain_plot = imread(fullfile([K_dir '/K' num2str(SelectK), '_3Dbrain.png']));
imwrite(cluster_brain_plot, fullfile(extract_results_dir, '3Dbrain.png'));

VectorNumber = imread(fullfile([K_dir '/K' num2str(SelectK), '_VectorNumber.png']));
imwrite(VectorNumber, fullfile(extract_results_dir, 'VectorNumber.png'));

ClusterPerformance = imread(fullfile([leida_res '/ClusterPerformance.png']));
imwrite(ClusterPerformance, fullfile(extract_results_dir, 'ClusterPerformance.png'));

K6_LinksCortex = imread(fullfile([K_dir '/K' num2str(SelectK), '_LinksCortex.png']));
imwrite(K6_LinksCortex, fullfile(extract_results_dir, 'K6_LinksCortex.png'));

ClusterPerformance = imread(fullfile([leida_res '/ClusterPerformance.png']));
imwrite(ClusterPerformance, fullfile(extract_results_dir, 'ClusterPerformance.png'));

Centroid_Pyramid_SideView = imread(fullfile([leida_res '/Centroid_Pyramid_SideView.png']));
imwrite(Centroid_Pyramid_SideView, fullfile(extract_results_dir, 'Centroid_Pyramid_SideView.png'));

% Data info 
Data_info = load([fullfile(leida_res,'/LEiDA_EigenVectors.mat')]).Data_info;
writetable(struct2table(Data_info), fullfile(extract_results_dir, 'data_info.csv'));

% Dwell time
LT = load([fullfile(leida_res,'/LEiDA_Stats_DwellTime.mat')]).LT;
LT = LT(:,SelectK - 1,1:SelectK); % Get dwell times for selected K
LT = reshape(LT, [length(LT), SelectK]);
writematrix(LT, fullfile(extract_results_dir, 'dwell_time.csv'));

% Fractional occupancy
P = load([fullfile(leida_res,'/LEiDA_Stats_FracOccup.mat')]).P;
P = P(:,SelectK - 1,1:SelectK); % Get dwell times for selected K
P = reshape(P, [length(P), SelectK]);
writematrix(P, fullfile(extract_results_dir, 'fractional_occupancy.csv'));

% Transition probaility
TMnorm = load(fullfile([K_dir '/LEiDA_Stats_TransitionMatrix.mat'])).TMnorm; %TM(s,alpha,beta) - alpha(columns) and beta(blocks) corresponds to departure and arrival states, respectively;

% Initialize a structure to store the matrices
matrices = struct();

% Extract each 12x6 slice from the 12x6x6 array and store in a structure
for i = 1:SelectK
    fieldName = sprintf('C%d', i);  % Create a field name dynamically (M1, M2, ..., M6)
    matrices.(fieldName) = TMnorm(:, :, i);
    writematrix( matrices.(fieldName) , fullfile(extract_results_dir, [fieldName, '_transitional_matrix.csv']));
end


