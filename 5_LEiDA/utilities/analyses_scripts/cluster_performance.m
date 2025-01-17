function cluster_performance(data_dir)
%
% Compute the Dunn's index, Calinski-Harabasz (CH) index and average 
% Silhouette coefficient for each value of K to assess clustering
% performance.
%
% INPUT:
% data_dir      directory where the results from K-means are saved
%
% OUTPUT:
% dunn_score    Dunn's index computed for each K
% avg_sil       average Silhouette coefficient computed for each K
% CH            CH index computed for each K
% .fig/.png     plot of the Dunn's index, average Silhouette
%               coefficient and CH index for each number of clusters
%
% Author: Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org
%         Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% Input example:
% data_dir = 'D:/LEiDA_Toolbox/LEiDA_results/';

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile(fullfile(data_dir, file_V1))  % CHANGED BY JL to make it more general
    load(fullfile(data_dir, file_V1), 'V1_all');
end
if isfile(fullfile(data_dir, file_cluster))
    load(fullfile(data_dir, file_cluster), 'Kmeans_results', 'rangeK');
end

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Dunn's index
% disp(' ')
% disp('Computing Dunn''s index:')
% distM = squareform(pdist(V1_all,'cosine'));
% dunn_score = zeros(length(rangeK),1);
% for i = 1:length(rangeK)
%     dunn_score(i) = dunns(rangeK(i), distM, Kmeans_results{i}.IDX);
%     disp(['- K = ' num2str(rangeK(i))])
% end
% [~, ind_maxdunn] = max(dunn_score);
% disp(['- Best clustering solution according to Dunn''s index: ' num2str(rangeK(ind_maxdunn)) ' clusters']);

% Randomly sample participants to calculate Dunn'a index, otherwise memory limit is exceeded
rows_per_participant = 193;
num_participants = size(V1_all, 1) / rows_per_participant;
num_sampled_participants = 100; % Adjust this number based on available memory
sampled_data = [];
sampled_IDX = cell(length(Kmeans_results), 1);

% Randomly sample participants
sampled_participants = randperm(num_participants, num_sampled_participants);

% Collect data and clustering indices for the sampled participants
for i = 1:num_sampled_participants
    participant = sampled_participants(i);
    start_idx = (participant - 1) * rows_per_participant + 1;
    end_idx = start_idx + rows_per_participant - 1;
    sampled_data = [sampled_data; V1_all(start_idx:end_idx, :)];
    
    % Sample the clustering indices for each K
    for k = 1:length(Kmeans_results)
        if isempty(sampled_IDX{k})
            sampled_IDX{k} = [];
        end
        sampled_IDX{k} = [sampled_IDX{k}, Kmeans_results{k}.IDX(start_idx:end_idx)];
    end
end

% Calculate pairwise cosine distances and convert to a square matrix
distM = squareform(pdist(sampled_data, 'cosine'));

% Initialize an array to store Dunn's index for each K
dunn_score = zeros(length(rangeK), 1);

% Loop over each K in the range
for i = 1:length(rangeK)
    % Calculate Dunn's index for the current number of clusters (rangeK(i))
    dunn_score(i) = dunns(rangeK(i), distM, sampled_IDX{i});
    % Display the current number of clusters being evaluated
    disp(['- K = ' num2str(rangeK(i))])
end

% Find the index of the maximum Dunn's index value
[~, ind_maxdunn] = max(dunn_score);

% Display the best number of clusters according to Dunn's index
disp(['- Best clustering solution according to Dunn''s index: ' num2str(rangeK(ind_maxdunn)) ' clusters']);

% Average Silhouette Coefficient
disp(' ')
disp('Computing average Silhouette coefficient:')
avg_sil = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_sil = evalclusters(V1_all,Kmeans_results{i}.IDX','Silhouette','Distance','cosine');
    avg_sil(i) = eva_sil.CriterionValues;
    clear eva_sil;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxsil] = max(avg_sil);
disp(['- Best clustering solution according to average Silhouette coefficient: ' num2str(rangeK(ind_maxsil)) ' clusters']);

% CH index
disp(' ')
disp('Computing CH index:')
CH = zeros(length(rangeK),1);
for i = 1:length(rangeK)
    eva_CH = evalclusters(V1_all,Kmeans_results{i}.IDX','CalinskiHarabasz');
    CH(i) = eva_CH.CriterionValues;
    clear eva_CH;
    disp(['- K = ' num2str(rangeK(i))])
end
[~, ind_maxCH] = max(CH);
disp(['- Best clustering solution according to CH index: ' num2str(rangeK(ind_maxCH)) ' clusters']);

% Saving results from clustering performance analysis
save_file = 'ClusterPerformance.mat';
save([data_dir '/' save_file],'CH','dunn_score','avg_sil');
disp(' ')
disp(['Clustering performance results saved successfully as ' save_file]);
disp(' ')

disp('Plotting clustering performance results:')
Fig = figure('Position', get(0, 'Screensize'));
x = 2:1:rangeK(end);
tiledlayout(3,1);
ax1 = nexttile;
plot(ax1,x,dunn_score,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax1,'Dunn''s index','Fontsize',12);
box off;

ax2 = nexttile;
plot(ax2,x,avg_sil,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax2,'Average Silhouette coefficient','Fontsize',12);
box off;

ax3 = nexttile;
plot(ax3,x,CH,'b','LineWidth',2);
xticks(2:1:rangeK(end));
ylabel(ax3,'CH index','Fontsize',12);
box off;
xlabel('Number of clusters','Fontsize',12);

saveas(Fig, fullfile(data_dir, 'ClusterPerformance.png'),'png');
saveas(Fig, fullfile(data_dir, 'ClusterPerformance.fig'),'fig');
disp('- Plot successfully saved as ClusterPerformance');

close all;
