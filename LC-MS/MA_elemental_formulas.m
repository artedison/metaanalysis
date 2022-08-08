%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Elemental Formulas for Meta-Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
x = readtable('rp_neg_analysisReady_sbys_wAnno_full.txt', 'Delimiter', '\t'); %load in merged and filtered annotation file and wide file containing ms2_id
y = readtable('MA_rp_neg_formula_identifications.txt', 'Delimiter', '\t'); %load in SIRIUS elemental formula summary file

%% Extract ms2 identifiers from each data matrix
id_x = string(x.ms2_id);
id_x_trunc = strtok(id_x, '_'); %remove collision energy information - leaving only a double

id_y = string(y.ms2_id); 

%% Find matching ms2_id between data matrices 
[Lia,Locb] = ismember(id_y,id_x_trunc); %output logical index and row index for the location of id_y found in id_x_trunc

Locb( ~any(Locb,2), : ) = [];  %remove empty rows in Locb

y_EF = y(Lia,:); %extract only the rows in y that were found in x

x_unqID = x.UniqueID; %create array containing the unique identifier for all features in X

metID_EF = x_unqID(Locb); %extract unique identifiers in X according to their Locb index

EF_table = [metID_EF y_EF]; %merge unique identifier with elemental formulas 
EF_table.Properties.VariableNames{1} = 'UniqueID';


%% Write to CSV 
writetable(EF_table,'SIRIUS_MA_rp_neg_mergedID.csv');

