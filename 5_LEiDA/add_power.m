% Load the .mat file
data = load('ParcelsMNI2mm.mat');

% Append additional data
parcels = niftiread('/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/power264MNI.nii');
labels = importdata('/Users/hannahcasey/Desktop/PhD/projects/STRADL_dynamic_rs_chronic_pain_depression/resources/tpl-MNI152NLin2009cAsym_atlas-power2011_dseg.txt');

charCellArray = repmat({''}, 264, 21);

% Initialize the 264x21 character array with spaces
charArray = char(zeros(264, 21) + ' ');

for i = 1:length(labels)
    % Convert each cell content to a character array and pad with spaces to ensure 21 characters
    paddedCharArray = pad(labels{i}, 21);
    
    % Assign the padded character array to the corresponding row in charArray
    charArray(i, :) = paddedCharArray;
end


data.Power.labels = charArray;
data.Power.volume = parcels;

% Save the updated structure back to the .mat file
save('ParcelsMNI2mm.mat', '-struct', 'data');