%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%  DIRECTORIES SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% OPTION 1: Manual
    % 
% 
% %             % Add project directory on YOUR local machine and adjust
% paths as necessary
%     
 [paths.matlab,paths.scripts,paths.thisWorkflow] = findCurrentFile(); 
 paths.project = ['/Users/gjg/Dropbox (Edison_Lab@UGA)/Projects/CIDC_worms/NMR_analysis/Sweet16_NMR']
 
%  %% Defining 
%  
%    % These is a standard directory structure
%         
        paths.results = [paths.project,'/','results'];
        paths.data = [paths.project,'/','data/D2O'];
            paths.raw = [paths.data,'/','raw/noesypr1d_batch1and2'];
            paths.processed_1_2 = [paths.data,'/','processed/noesypr1d_batch1and2'];             
               paths.processed_3_4 = [paths.data,'/','processed/noesypr1d_batch3and4'];    
                  paths.processed_5_6 = [paths.data,'/','processed/noesypr1d_batch5and6'];    
            paths.metadata = [paths.project,'/','metadata'];     
        paths.design_file = [paths.project,'/','design_file'];  
%         
%         addpath(genpath(paths.project))
%  
%         %%
%        
       
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    DATA IMPORT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import DATA 

  
  % Importing batch1 and 2
cd (paths.processed_1_2)
cd ft_edit

spectra_1_2 = loadallft_GG

  % Importing batch3 and 4
cd (paths.processed_3_4)
cd ft_edit

spectra_3_4 = loadallft_GG

  % Importing batch5 and 6
cd (paths.processed_5_6)
cd ft_edit

spectra_5_6 = loadallft_GG
clearvars spectra 
%% Merging all datasets into one variable 

spectra_all = [spectra_1_2...
    spectra_3_4...
    spectra_5_6]

clearvars spectra_5_6 spectra_1_2 spectra_3_4 


%% Import Metadata
cd (paths.metadata)
T = readtable('master_nmr_run - Sheet1.csv') % this contains all the data including D2O and CHCl3 

Td = T(contains(T.Solvent,"D2O"),:) %selecting only D2O



%%        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    DATA IMPORT    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %% Referencing/Calibration (e.g. to DSS)

 %%
    % Manual peak selection (set for DSS by default)
    
 
referenced_spectra= ref_spectra_manual (spectra_all, [-0.2, 0.2], 0.01,0.00)

        
 
      %% Convert data to spectral matrix ("X matrix")
    

    % Make the X matrix and ppm vector
    
        [wrk_data.X,wrk_data.ppm,~]=Setup1D(referenced_spectra);    
        
    % Plot the matrix
       
        figure, plotr(wrk_data.ppm,wrk_data.X)
         
            title('Referenced Spectra')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            whichLine()
            %fig = gca;
            %saveas(fig,fig.Title.String,'png')
          
       % run this after the plot is up  to select spectra that will need to be rephased      
% [wrk_data.rephase1] = selectOutliers()






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%         ALIGNMENT        %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Usage:
%% Pre-align spectra and build prealigns input structure (gives user more flexibility with alignment params)
    currentppm = wrk_data.ppm;
    matrix =wrk_data.X;

    % Suppress figures:
        set(0,'DefaultFigureVisible','off')
    % Plot the unaligned data for comparison:
             wrk_data.alignment.preAligns(1).Data = matrix;
             wrk_data.alignment.preAligns(1).Titles = 'Unaligned';

                 % Try the Pearson PAFFT
        % Plot it:
            wrk_data.alignment.preAligns(2).Data = guide_align1D(matrix,currentppm,'correlation','CCOW');
            wrk_data.alignment.preAligns(2).Titles = 'XRALg - Guide Align - CCOW (correlation)';
 
   % Try the Spearman PAFFT     
        % Plot it:
              wrk_data.alignment.preAligns(3).Data = star_align1D(matrix,currentppm,'mean','RAFFT' );
              wrk_data.alignment.preAligns(3).Titles = 'XRALg - Guide Align - RAFFT ()';
            
 
        set(0,'DefaultFigureVisible','on') % turn figures back on
        clear('currentppm','matrix')
        
%% Save a copy of the workspace 

    save('phased_merged_aligned_3Aug2020.mat')
        
        
%% Compare alignments, choose region-by-region, and gather consensus 
  cd (paths.results)
  cd sweet_16_NMR_d2o_3Aug2020
    wrk_data.alignment.outdir = 'output_multiAlign';
    wrk_data.alignment.lastFig = []; % comment out if re-running from previous figure
    wrk_data.alignment.patches = {}; % comment out if re-running from previous figure
    
    [wrk_data.XRBA,wrk_data.alignment.multiAlign_out] = multiAlign_interactive_draft_2(wrk_data.X,...
                                                           wrk_data.ppm,...
                                                           wrk_data.alignment.outdir,...
                                                           wrk_data.alignment.preAligns,...
                                                           wrk_data.alignment.lastFig,...
                                                           wrk_data.alignment.patches);
    paths.alignments = cd();                                                 
    
%% Clear alignment vars, save
     save('post_alignment_03Aug2020.mat') % after running, comment this out so we
                                 % don't accidentally overwrite it at some point
    % load('post_alignment.mat') % ^ (same here). 
    
    %% 
    
%% Picking regions to remove 

    % Select the region
        figure, plot(wrk_data.ppm,wrk_data.XRBA)
            set(gca,'XDir','reverse')
            title('Select the Region to Remove (i.e. Solvent peak)')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
        [wrk_data.solvent_region] = selectROIsFromFigure();close(gcf)
        
    % Skip selection
     
      
%%
wrk_data.XR=wrk_data.XRBA;% storing as a new variable for data recovery

   % % Remove Selected regions , H2O etc
     for i= 1:length(wrk_data.solvent_region)
         
        wrk_data.XR = remove_region_GG(wrk_data.XR,wrk_data.ppm,wrk_data.solvent_region(1,i),wrk_data.solvent_region(2,i));
       
     end
     
     
     [wrk_data.XR,wrk_data.ppmR]=remove_ends(wrk_data.XR,wrk_data.ppm,0.5,9.66)
     
     
        figure, plot(wrk_data.ppmR,wrk_data.XR)
            set(gca,'XDir','reverse')
            title('Trimmed Spectra (no water signal)')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            whichLine()

            
          
 %% Removing blanks and getting a matrix of just blanks

wrk_data.Xn_blks = wrk_data.XR(find(~contains( Td.genotype,'d20')),:) %making a blank matrix

Tdata = Td(find(~contains( Td.genotype,'d20')),:) %making a metadata table for just data (blanks contain d20 in the name)

wrk_data.X_blanks = wrk_data.XR(find(contains( Td.genotype,'d20')),:) % making table for just blanks 

figure,
plotr(wrk_data.ppmR,wrk_data.Xn_blks)

figure,
plotr(wrk_data.ppmR,wrk_data.XR(find(contains( Td.genotype,'d20')),:))
 %% Baseline Correction

%get the baseline factors A and S
    [wrk_data.baselineCorr.A,wrk_data.baselineCorr.S] = showBaseline(wrk_data.Xn_blks(11,:))
   
 %%
 
 % apply and visualize correction
    wrk_data.baselineCorr.A =   1e+07
    wrk_data.baselineCorr.S = 3000
    showBaseline (max(wrk_data.Xn_blks), wrk_data.ppmR, wrk_data.baselineCorr.A, wrk_data.baselineCorr.S);
    
   %%
   %apply to the whole data set
    wrk_data.XRB = CorrectBl(wrk_data.Xn_blks,      wrk_data.baselineCorr.A,      wrk_data.baselineCorr.S);
   
    
    %%
    
  
    figure, 
    plotr(wrk_data.ppmR,wrk_data.XRB)
           
            title('Baseline-Corrected Spectra')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
    %    [wrk_data.outlies] = selectOutliers()
            

%%

    
    %% colors by a table field
    
      figure,   
    plotr_colorby(wrk_data.ppmR,wrk_data.XRB,Td.genotype)
           
            title('Aligned spectra')
            xlabel('Chemical Shift (ppm)')
            ylabel('Signal Intensity')
            
            
            
    %% Setting up a matrix for stacked spectra in group order
    
[test_sort, ind_test] = sort(Td.genotype (~contains(Td.genotype,'d20'),:))

[C,ia,ic]= unique(test_sort)

sortdata = wrk_data.XRB(ind_test,:);

%%
stackSpectra_2(sortdata,wrk_data.ppmR,0.011,0.001,'white','lines',ic,[])

clearvars C ia ic sortdata test_sort ind_test
%% 


    %% Optimize Peak Picking (threshold; for representative spectrum)


        matrix = wrk_data.XRB;
        ppm = wrk_data.ppmR;

             optimize_Peakpick1D(matrix,ppm,'var',0.05:0.05:0.5,'Complex');   

             doToAllOpenFigs('set(gca,''xlim'',[0.7,1.7])')    
               

        clear('matrix','ppm')    
        
        

    %% Do the actual Peak Picking (for representative spectrum)

        matrix = wrk_data.XRB;
        ppm = wrk_data.ppmR;

             peaks = struct();
        
             [peaks.ints, peaks.shifts]= Peakpick1D(matrix ,ppm,'max',0.2,'Complex');

  clearvars matrix ppm 
             %%
     figure,
     hold
     
     plotr(wrk_data.ppmR,wrk_data.XRB, 'r','LineWidth',1.5)
     plotr(peaks(1).shifts,peaks(1).ints, 'co')

         
%% Automatic Binning (Bucketing)


%% Generate buckets using a range of both params

        matrix = wrk_data.XRB;
        ppm = wrk_data.ppmR;
    
    sb = 0.002:0.002:0.008;
    sl = 0.3:0.1:0.6;
  
    [optOB_out] = optimize_optBucket(matrix,ppm,sb,sl);

     
    %% Filter out the bins with no peaks
    
       [optOB_out] = filterBuckets_Peaks_opt(ppm,optOB_out, peaks);        
 
    %% Plot the results of optOB

         [optOB_out] =  plotOptBucket_optResult(matrix,ppm,optOB_out,[3.6584    4.0], [min(matrix(:)) max( matrix(:,(ppm>3.6584  & ppm<4.0)), [],'all' ) ]);


             %% Peak pick every spectra independently, so that each spectrum contributes to the peak shape between the boundaries
             matrix = wrk_data.XRB;
        ppm = wrk_data.ppmR;
        peakthresh = 0.15
        mode = 'Complex'
        
             
              [itpeaks] = Peakpick1D_per_spectra(matrix,ppm,peakthresh,mode)


        
 %%  Expand buckets to each of the bins boundaries
 
 
    [optOB_out] = expandBucketBounds (optOB_out, matrix, ppm, 'plotResult');
    
    
            clearvars matrix ppm     
%%  Manual Refinement of the boundaries
    %% [optOB_out] = refineBuckets(matrix,ppm,optOB_out,5);    
    
        matrix = wrk_data.XRB;
        ppm = wrk_data.ppmR;
        
     [optOB_out] = refineBuckets_GG2(matrix,ppm,optOB_out,itpeaks,peaks, 'expandedBuckets');    

     %%
       save('post_buckets_3Aug2020.mat')
%%
%% This calculates the number of peaks in each bin in the original Peakpick1D_per_spectra data
        % Calculates the max peak within each bin and its chemical shift 
        % gap fills both ppm and intensities of peaks that were no present
        % in the bin or not picked (not detected/below baseline)

matrix = wrk_data.XRB;
ppm = wrk_data.ppmR; 
 buckets = optOB_out.refinedBuckets.refinedBuckets;
 ppeak_struct = itpeaks;
    
        
[itpeaks] = peaks_per_bin (matrix, ppm,buckets, ppeak_struct);


%% Align metric

matrix = wrk_data.XRB;
ppm = wrk_data.ppmR; 
 buckets = itpeaks.max_iteration.sorted_bucket_list;
 ppeak_struct = itpeaks;
 metadata = Tdata;
 
itpeaks = align_metric ( matrix, ppm, ppeak_struct, buckets, metadata);

%% getting blank matrix peak picked

blank_matrix = wrk_data.X_blanks;
ppm = wrk_data.ppmR; 
[peaks.blank_ints, peaks.blank_shifts]= Peakpick1D(blank_matrix ,ppm,'max',0.9,'Complex');

%% scoring blank peaks 

spectra_struct = wrk_data.XRBA;
ppick_output = peaks;
per_spectra_ppick_output = itpeaks;
sorted_bucket_list = itpeaks.max_iteration.sorted_bucket_list;

 [itpeaks] = blank_feature_score (spectra_struct, blank_matrix, ppm, ppick_output,  per_spectra_ppick_output, sorted_bucket_list)


%% storing matrices
matrix = itpeaks.max_iteration.gap_f_ints ;
ppm = itpeaks.max_iteration.gap_f_median_ppm	; 
metadata_1 = Tdata.wormgrowth_sample_name;
metadata_2 = Tdata.run_order;
prefix = 'ph';

ppm_string = join(...
    horzcat(...
    repmat({'ppm'},length(ppm),1),...
    replace(string(ppm), '.','_')...
    ),...
    '_',2);

Varnames = [  {'blank_peak'};  {'blank_intensity_flag'}; {'align_score'}; join(...
    horzcat(...
                string(metadata_1),...
                repmat( {'NMR'}, length(metadata_1),1),...
                repmat( {'D2O'}, length(metadata_1),1),...
                string(metadata_2)),...
                '_') ];

table_matrix = array2table( ...
    [ppm_string,...
    [ itpeaks.blanks.blanks_flag,  itpeaks.blanks.blank_ints_flag', itpeaks.align_scale.metric', matrix']],...
    'VariableNames' ,['ppm';Varnames]); 

writetable (table_matrix, join([string(prefix), '_NMR_D2O_14Aug2020.txt'],''), 'WriteRowNames', true, 'Delimiter', '\t')

clearvars matrix ppm metadata_1 metadata_2 ppm_string  table_matrix
%%
matrix = itpeaks.max_iteration.gap_f_ints ;
ppm = itpeaks.max_iteration.gap_f_median_ppm	; 
metadata_1 = Tdata.wormgrowth_sample_name;
metadata_2 = Tdata.run_order;
prefix = 'ph';

   sampleID = array2table(join(...
    horzcat(...
                string(metadata_1),...
                repmat( {'NMR'}, length(metadata_1),1),...
                repmat( {'D2O'}, length(metadata_1),1),...
                string(metadata_2)),...
                '_'),...
                'VariableNames' , {'sampleID'})
    
    Tdwrite = horzcat( sampleID, Tdata);
    
    writetable( Tdwrite, 'Design_file_NMR_D2O_14Aug2020.txt', 'Delimiter', '\t')

    clearvars matrix ppm metadata_1 metadata_2  
    
    
