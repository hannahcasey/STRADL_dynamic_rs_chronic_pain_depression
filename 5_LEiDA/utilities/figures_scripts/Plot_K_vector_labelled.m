function Plot_K_vector_labelled(data_dir,save_dir,selectedK,centroid,parcellation)
%
% Plot the selected centroid in vector format with the labels of each area.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% centroid      centroid defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot centtroid as a barplot and label each area
%               according to the applied parcellation
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile(fullfile(data_dir, file_cluster))% changed by JL to make it more general
    load(fullfile(data_dir, file_cluster), 'Kmeans_results', 'rangeK');
end

% Matrix (selectedK*n_areas) where each row represents one PL state
V = Kmeans_results{rangeK == selectedK}.C;
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));
clear Kmeans_results

% Get the number of areas based on centroids size
n_areas = size(V,1);

% Get labels from parcellation used
% % commented by JL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %switch parcellation % commented by JL
% switch parcellation.name
%     case 'AAL116'
%         %load('ParcelsMNI2mm','label116'); 
%        %labels = label116([1:n_areas],:); % commented by JL
%        labels = label116(parcellation.areas,:); % added by JL
%         %clear label116
%     case 'AAL120'
%         load('ParcelsMNI2mm','label120');
%         %labels = label120([1:n_areas],:); % commented by JL
%         labels = label120(parcellation.areas,:); % added by JL
%         clear label120    
%     case 'brainnetome_cortical_atlas_yeo_dimensions'
%         load('ParcelsMNI2mm','brainnetome_cortical_labels_210');
%         %labels = brainnetome_cortical_labels_210([1:n_areas],:); % commented by JL
%         labels = brainnetome_cortical_labels_210(parcellation.areas,:); % added by JL
%         clear brainnetome_cortical_labels_210
% end
% % commented by JL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % added by JL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels_struct  = load('ParcelsMNI2mm', parcellation.name); % added by JL
labels = labels_struct.(parcellation.name).labels; % added by JL
labels = labels(parcellation.areas,:); % added by JL
clear labels_struct

% Reorder the position of each parcel
%Order = [1:2:n_areas n_areas:-2:2]; % commented by LK 18/02/2023
Order = 1:n_areas; % added by LK

disp(' ');
disp(['Plotting the vector format with labels for the centroid ' num2str(centroid) ':'])
Fig = figure('Position', get(0, 'Screensize')); 
subplot(1,3,2)
Vo = V(Order,centroid);
hold on
barh(Vo.*(Vo < 0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
barh(Vo.*(Vo >= 0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
ylim([0 n_areas+1])
xlim([-1 1])
set(gca,'YTick',1:n_areas)
set(gca,'YTickLabel',labels(Order,:),'Fontsize',6)
ax = gca;
ax.XAxis.FontSize = 10;
grid on
title(['PL state ' num2str(c)],'Fontsize',12)

saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_VectorLabels.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_VectorLabels.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) 'C' num2str(centroid) '_VectorLabels']);

close all;