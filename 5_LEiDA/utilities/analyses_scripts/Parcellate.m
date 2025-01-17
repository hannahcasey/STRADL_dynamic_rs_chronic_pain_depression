function Parcellate(data_dir,save_dir,parcellation)
%
% Parcellate the neuroimaging data according to a given parcellation.
% Accepts only .nii files.
%
% INPUT:
% data_dir      directory where the .1D files with the ABIDE data are
%               stored
% save_dir      directory to save the new data
% parcellation  parcellation atlas to be used to parcel the data
% parcellation.n_areas       number of areas of the parcellation to consider
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

% INPUT EXAMPLES:
% data_dir = 'D:/LEiDA_Toolbox/Outputs/dparsf/nofilt_noglobal/';
% save_dir = 'D:/LEiDA_Toolbox/ABIDE_dparsf_AAL120/';
% parcellation = 'AAL120'; (AAL116, AAL120, dbs80, glasser360)
% parcellation.n_areas = 94;

%%%%% commented by JL
% % Get number of files in folder
% % aux_data = [dir(fullfile([data_dir '*.nii'])); dir(fullfile([data_dir '*.nii.gz']))];
% aux_data = dir([data_dir '*.nii.gz']);
% num_subjs = numel(aux_data);
% 
% % Unzip data (comment if already unzipped)
% gunzip([data_dir '*.nii.gz']);
% num_subjs = numel(aux_data);

%%%%% commented by JL

%%%% added by JL
% check if there are unzipped nii in folder
aux_data = dir(fullfile(data_dir, '*.nii'));

% if not, check for zipped files and unzip them.
if isempty(aux_data)
    aux_data = dir(fullfile(data_dir, '*.nii.gz'));
    gunzip(aux_data); % added by JL
end

% check how many files there are
num_subjs = length(aux_data);
%%% added by JL





% Obtain parcellation atlas from ParcelsMNI2mm (Dr. Joana Cabral)
% V_Parcel = struct2array(load('ParcelsMNI2mm',['V_' parcellation])); %commented by JL
V_Parcel = struct2array(load('ParcelsMNI2mm',parcellation.name)).volume; % added by JL
sz = size(V_Parcel); % size of the parcellation

for s = 1:num_subjs
    
    file = aux_data(s).name;
    [~, baseFileName,~] = fileparts(file);
    
    if size(file,1)
        % disp(['Parcellating data file ' baseFileName ' using parcellation ' parcellation]); %commented by JL
        disp(['Parcellating data file ' baseFileName ' using parcellation ' parcellation.name]); % added by JL
        % Read the nii file
        %fMRI_MNI = niftiread([data_dir file]);  %commented by JL
        fMRI_MNI = niftiread(fullfile(data_dir, file));% added by JL
        T = size(fMRI_MNI,4); % number of volumes
        
        % Check if nii files correspond to MNI2mm
        if size(fMRI_MNI,1) ~= sz(1) || size(fMRI_MNI,2) ~= sz(2) || size(fMRI_MNI,3) ~= sz(3)
            % Files will be resized in order to be accomodated to the MNI2mm template
            fMRI_MNI2mm = zeros(sz(1), sz(2), sz(3), T);
            for t = 1:T
                fMRI_MNI2mm(:,:,:,t) = imresize3(fMRI_MNI(:,:,:,t),sz);
            end
        end
                
        fMRI_parcel = zeros(parcellation.n_areas,T);
        
        for n = 1:parcellation.n_areas
            ind_voxels = find(V_Parcel == n);
            
            for v = 1:numel(ind_voxels)
                [I1,I2,I3] = ind2sub(sz,ind_voxels(v));
                if ~isnan(fMRI_MNI2mm(I1,I2,I3,1))
                    fMRI_parcel(n,:) = fMRI_parcel(n,:) + squeeze(fMRI_MNI2mm(I1,I2,I3,:))';
                end
            end

            fMRI_parcel(n,:) = fMRI_parcel(n,:)/numel(ind_voxels);
            fMRI_parcel(n,:) = detrend(fMRI_parcel(n,:) - mean(fMRI_parcel(n,:)));
            fMRI_parcel(n,:) = fMRI_parcel(n,:)/std(fMRI_parcel(n,:));
        end
        
       % save([save_dir baseFileName '_' parcellation], 'fMRI_parcel') % commented by JL
        save(fullfile(save_dir, [baseFileName '_' parcellation.name]), 'fMRI_parcel') % added by JL
    end
end

% Unzipping the file if needed
% tmpDir = tempname;
% mkdir(tmpDir);
% if strcmp(file(end-2:end), '.gz')
%     [~, baseFileName,~] = fileparts(file(1:end-3)); % remove gz
%     try
%         tmpFile = gunzip([data_dir file], tmpDir);
%         % Read the nii file
%         fMRI_MNI = niftiread(char(tmpFile));
%     catch me
%     end
% else
%     [~, baseFileName,~] = fileparts(file);
%     % Read the nii file
%     fMRI_MNI = niftiread([data_dir file]);
%     rmdir(tmpDir);
%     % Read the nii file
% end