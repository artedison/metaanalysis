
%% Load data files 
x = readtable(, 'Delimiter', ',',  'ReadVariableNames',false); %matrix containing file names and strain information
y = readtable('BFF_corrected_rp_neg_with_rt_mz.txt', 'Delimiter', '\t'); %matrix containing peak intensities and file names


%% Extract only PD1074 samples from Data Matrix 

y_names = y.Properties.VariableNames; % pull out header names from data matrix 'y' 
strains = string(x.strain); % pull out strain names from data matrix 'x'
strains_pd1074 = contains(strains,'PD1074'); %logical vector identifying strain 'PD1074'

pd1074_filenames = unique(x.mzml_file_name(strains_pd1074)); %vector containing file names associated with strain 'PD1074'

y_pd1074_logical = contains(y_names,pd1074_filenames); %logical vector identifying which column names in 'y' match with 'PD1074' filenames

y_pd1074_subset = y(:,y_pd1074_logical); %subset of data matrix 'y' with only PD1074 strains
%% Identify stable features (features in 100% of individual PD1074 Control Strain)
 
z = table2array(y_pd1074_subset); %convert 'y_pd1074_subset' from a table to an array
z_zeros_logical = z>0; %identify features with '0' as intensity value
z_zeros_logical_sum = sum(z_zeros_logical,2); %sum rows 

z_100pds_logical = z_zeros_logical_sum==size(z,2); %identify features that have a value > 0 in all PD1074 strains  (100%)

y_100pd_subset = y(z_100pds_logical,:); %Subset original data matrix 'y' with only features found in 100% of PD1074 control strain


%% remove features based on retention time 
rt_min = 0.65;
rt_max = 9.0;

rt_log = (y_100pd_subset.rt < rt_min) | (y_100pd_subset.rt > rt_max);
sum(rt_log)

y_100pd_subset_rt = y_100pd_subset(rt_log==0,:);

%% Write to table and export
writetable(y_100pd_subset_rt,'BFF_corrected_rp_neg_with_rt_mz_solv_filt.txt');
