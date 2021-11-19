% some preliminaries for MATLAB path
clear 
close all
clc

% initialize fieldtrip
path_ft   = 'D:\Dokumentumok\MATLAB\fieldtrip\fieldtrip-20211015'; %path for the fieldtrip folder
path_data = 'D:\Dokumentumok\MATLAB\EEG\Data'; %path for the data

% adding fieldtrip to the matlab path + the relevant subdirectories
addpath(path_ft);
ft_defaults;

%% DECLARATION of constants that are important for future analysis
% IDs for all subjects, then divided into experimental and control group
% for possible need in the future
subjects_all = {'a564-j168' 'b574-r862' 'c472-t464' 'e636-d750' 'e814-a939' 'f923-x353' 'g197-r432' 'g627-x634' 'g977-e686' 'm927-k134' 'n623-r167' 'n949-f538' 'o565-o662' 'p293-a263' 'r437-k431' 'r684-b222' 'r967-g952' 's422-v155' 'u184-p978' 'u338-f143' 'v146-f150' 'v511-b144' 'a179-m176' 'b277-c676' 'b436-e194' 'k692-n675' 'l434-n913' 'n476-v576' 'n535-e890' 'n683-j492' 'n713-j980' 'o391-h387' 'o522-h323' 'p359-s679' 'p665-c271' 'q511-e393' 'q543-e754' 'q688-x553' 's493-s280' 's979-u715' 't473-l690'};
subjects_exp = {'a564-j168' 'b574-r862' 'c472-t464' 'e636-d750' 'e814-a939' 'f923-x353' 'g197-r432' 'g627-x634' 'g977-e686' 'm927-k134' 'n623-r167' 'n949-f538' 'o565-o662' 'p293-a263' 'r437-k431' 'r684-b222' 'r967-g952' 's422-v155' 'u184-p978' 'u338-f143' 'v146-f150' 'v511-b144' };
subjects_cont = {'a179-m176' 'b277-c676' 'b436-e194' 'k692-n675' 'l434-n913' 'n476-v576' 'n535-e890' 'n683-j492' 'n713-j980' 'o391-h387' 'o522-h323' 'p359-s679' 'p665-c271' 'q511-e393' 'q543-e754' 'q688-x553' 's493-s280' 's979-u715' 't473-l690'};

%the variable conditions contains the name of the 20 conditions we had
conditions = {'trans_stand' 'trans_dev' 'melody_stand' 'melody_dev' 'mistuning_stand' 'mistuning_dev' 'rhythm_short_stand' 'rhythm_short_dev' 'rhythm_long_stand' 'rhythm_long_dev' 'timing_mid_stand' 'timing_mid_dev' 'timing_end_stand' 'timing_end_dev' 'timbre_mid_stand' 'timbre_mid_dev' 'timbre_end_stand' 'rhythm_stand' 'rhythm_dev' 'timing_stand' 'timing_dev'};

%important triggers matched with the conditions (20 set of triggers, most
%of them contains only one trigger)
triggervalues = {'S  1' 'S  5' 'S  3' 'S  6' 'S  3' 'S  7' 'S  2' 'S 10' 'S  3' 'S101' 'S  3' 'S 11' 'S  4' 'S 12' 'S  3' 'S  8' 'S  4' {'S  2' 'S  3'} {'S 10' 'S101'} {'S  3' 'S  4'} {'S 11' 'S 12'}};

%% DATA EPOCHING for all subjects for all conditions
tic %starting timer just for our information
for j=1:length(conditions) %going through all 20 conditions
    for i=1:length(subjects_all) %going through every subject (in every condition)
        %saving the path and name of the EEG files that are needed to be
        %into a varible named 'dataset'
        %the data were exported from Brain Vision Analyzer(BVA), all
        %preprocessing is done in that software
        dataset  = fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', conditions{j} ,'.dat'));
        
        % in this step the epochs-of-interest are defined
        cfg           = [];
        cfg.dataset   = dataset;
        cfg.trialfun  = 'ft_trialfun_general'; 
        cfg.trialdef.eventtype = 'Stimulus';
        cfg.trialdef.eventvalue = triggervalues{j}; %from triggervalues variable, the index is matched as I defined the variable to do so
        cfg.trialdef.prestim    = 0.1; %0.1 secs prestimulus
        cfg.trialdef.poststim   = 0.7; %0.7 secs poststimulus
        cfg           = ft_definetrial(cfg);
        trl           = cfg.trl; % this field contains the definition of the epochs of interest 

        cfg.outputfile = fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', conditions{j} ,'_preprocessed')); %giving the name for the outputfile
        ft_preprocessing(cfg); %no other parameters defined, as preprocessing(eg. ICA, baseline-correction, artifact rejection) was done in BVA
    end
end
toc % Elapsed time is .... seconds.

%% COMPUTATION OF EVENT-RELATED POTENTIALS AND OF THE GROUP AVERAGE
% only for one condition: MISTUNING
%for mistuing standard individual average for pariticipants
tic
for i=1:length(subjects_all) %going through all subjects
    %decalre a variable for the complicated filename
    datafile = fullfile(path_data, strcat('music_mfeat_', subjects_all{i} , '_3_', 'mistuning_stand' ,'_preprocessed.mat'));
    load(datafile);
        
    % average across trials for all channels
    cfg         = [];
    erp         = ft_timelockanalysis(cfg,data);

    %saving the average ERPs into individual files
    save(fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', 'mistuning_stand' , '_erp.mat')), 'erp');
end
toc

% for mistuning deviant individual average for pariticipants
tic
for i=1:length(subjects_all)
    datafile = fullfile(path_data, strcat('music_mfeat_', subjects_all{i} , '_3_', 'mistuning_dev' ,'_preprocessed.mat'));
    load(datafile);
        
    % average across trials
    cfg         = [];
    erp         = ft_timelockanalysis(cfg,data);
    save(fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', 'mistuning_dev' , '_erp.mat')), 'erp');
end
toc

%% load previously calculated single subject ERPs in a single cell array - MISTUNING
% cell array: datatype
% -indexed data
% eg 0x0
allsubj_mistuning_stand = cell(1,length(subjects_all)); %creating cell array 'allsubj_mistuning_stand' and fill in with the for loop
for i=1:length(subjects_all)

    datafile     = fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', 'mistuning_stand' , '_erp.mat'));
    load(datafile);
    allsubj_mistuning_stand{i} = erp; %erp is a variable in the loaded file
end

allsubj_mistuning_dev = cell(1,length(subjects_all)); %creating cell array and fill in with the for loop
for i=1:length(subjects_all)

    datafile     = fullfile(path_data, strcat('music_mfeat_', subjects_all{i}, '_3_', 'mistuning_dev' , '_erp.mat'));
    load(datafile);
    allsubj_mistuning_dev{i} = erp;
end

% save the two cell arrays into a file, so that we do not need to run all
% the analysis from the start
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning_ind' , '.mat')), 'allsubj_mistuning_stand', 'allsubj_mistuning_dev');

%% calculate GROUP AVERAGE for MISTUNING standard and deviant

% load individual subject data that was previously saved 
load(fullfile(path_data, strcat('music_mfeat_', 'mistuning_ind' , '.mat')));

% calculate grand average for both condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grandavg_mistuning_stand  = ft_timelockgrandaverage(cfg, allsubj_mistuning_stand{:});
grandavg_mistuning_dev   = ft_timelockgrandaverage(cfg, allsubj_mistuning_dev{:});
% "{:}" means to use data from all elements of the variable

%% VISUALIZATION OF THE GROUP AVERAGE ERP

% % plot a sensor layout with the ERP time course at each electrode
cfg = [];
cfg.showlabels  = 'yes';
cfg.showscale = 'yes';
cfg.layout      = 'acticap-64ch-standard2.mat'; %this is our EEG systems's layout
figure;
ft_multiplotER(cfg, grandavg_mistuning_stand, grandavg_mistuning_dev);

% plot one chosen channel
cfg = [];
cfg.showlabels  = 'yes';
cfg.showlegend = 'yes';
cfg.channel = 'Fz';
cfg.linewidth  = 2;
figure;
ft_singleplotER(cfg, grandavg_mistuning_stand, grandavg_mistuning_dev);

% topoplot standard at 240 ms
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.xlim = [0.24 0.24]; % to define the timepoint
cfg.zlim = [-3 0.5]; %to have a uniform colorbar scale
cfg.parameter = 'avg'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,grandavg_mistuning_stand); colorbar


% topoplot deviant at 240 ms
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.xlim = [0.24 0.24];
cfg.zlim = [-3 0.5];
cfg.parameter = 'avg'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,grandavg_mistuning_dev); colorbar

%% difference of standard and deviant
% apply ft_math
% subtract avgCondition1 from avgCondition2

cfg = [];
cfg.operation    = 'subtract';
cfg.parameter    = 'avg'; %ERP
dev_minus_stand = ft_math(cfg, grandavg_mistuning_dev, grandavg_mistuning_stand);


% singleplot difference one chosen channel
cfg = [];
cfg.showlabels  = 'yes';
cfg.showlegend = 'yes';
cfg.channel = 'Fz';
cfg.linewidth  = 2;
figure;
ft_singleplotER(cfg, dev_minus_stand);

% multiplot of the difference
cfg = [];
cfg.showlabels  = 'yes';
cfg.showscale = 'yes';
cfg.layout      = 'acticap-64ch-standard2.mat'; %this is our EEG systems's layout
figure;
ft_multiplotER(cfg, dev_minus_stand);

% topoplot the difference for the whole period
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.xlim = [0.24 0.24];
cfg.zlim = [-3 0.5];
cfg.parameter = 'avg'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,dev_minus_stand); colorbar


%% T-TEST with FieldTrip function for ONE CHANNEL
% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'Fz';
cfg.latency     = [0.2 0.3];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';

Nsub = length(subjects_all);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubj_mistuning_stand{:}, allsubj_mistuning_dev{:});   % don't forget the {:}
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_t-test_Fzchan', '.mat')), 'stat');

%% PARAMETRIC T-TEST with MATLAB function for ALL CHANNELS without correction
%loop over channels
% time = [0.2 0.3];
% timesel_mistuning_stand = find(grandavg_mistuning_stand.time >= time(1) & grandavg_mistuning_stand.time <= time(2));
% timesel_mistuning_dev  = find(grandavg_mistuning_dev.time >= time(1) & grandavg_mistuning_dev.time <= time(2));
% clear h p
% 
% mistuning_stand_min_mistuning_dev = zeros(1,length(subjects_all));
% 
% for iChan = 1:65
%     for isub = 1:length(subjects_all)
%         mistuning_stand_min_mistuning_dev(isub) = ...
%             mean(allsubj_mistuning_stand{isub}.avg(iChan,timesel_mistuning_stand)) - ...
%             mean(allsubj_mistuning_dev{isub}.avg(iChan,timesel_mistuning_dev));
%     end
% 
%     [h(iChan), p(iChan)] = ttest(mistuning_stand_min_mistuning_dev, 0, 0.05 ); % test each channel separately
% end

%% PARAMTERIC T-TESTS WITH CORRECTION
%%Bonferroni correction with FieldTrip function

cfg = [];
cfg.channel     = 'all'; %now all channels
cfg.latency     = [0.2 0.3]; %predefined, by looking at the ERP-plots; MMN -> ~150-250 ms
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'bonferroni';

Nsub = length(subjects_all);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubj_mistuning_stand{:}, allsubj_mistuning_dev{:});
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_t-test_allchan_bonferroni', '.mat')), 'stat');

% plot corrected significant channels
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'acticap-64ch-standard2.mat';
cfg.highlight = 'on';
cfg.highlightcolor = 'r';
cfg.highlightchannel = find(stat.prob<(0.05/65));
cfg.comment   = 'no';
figure;
ft_topoplotER(cfg, grandavg_mistuning_stand)
title('Significant channels with multiple comparison correction')

%%corrected multiple comparisons with MATLAB function for ALL CHANNELS
%loop over channels
% time = [0.2 0.3];
% timesel_mistuning_stand = find(grandavg_mistuning_stand.time >= time(1) & grandavg_mistuning_stand.time <= time(2));
% timesel_mistuning_dev  = find(grandavg_mistuning_dev.time >= time(1) & grandavg_mistuning_dev.time <= time(2));
% clear h p
% 
% mistuning_stand_min_mistuning_dev = zeros(1,length(subjects_all));
% 
% for iChan = 1:65
%     for isub = 1:length(subjects_all)
%         mistuning_stand_min_mistuning_dev(isub) = ...
%             mean(allsubj_mistuning_stand{isub}.avg(iChan,timesel_mistuning_stand)) - ...
%             mean(allsubj_mistuning_dev{isub}.avg(iChan,timesel_mistuning_dev));
%     end
% 
%     [h(iChan), p(iChan)] = ttest(mistuning_stand_min_mistuning_dev, 0, 0.05/65 ); % test each channel separately
% end

% plot corrected "significant" channels
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = 'acticap-64ch-standard2.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(h);
% cfg.comment   = 'no';
% figure;
% ft_topoplotER(cfg, grandavg_mistuning_stand)
% title('Significant channels with multiple comparison correction')

%% NON-PARAMETRIC TESTS - more advanced
% montecarlo
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [0.2 0.3];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 100000;  % there are 41 subjects, so 2^41 =2,199,023,255,552 possible permutations, don't need all

Nsub = length(subjects_all);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubj_mistuning_stand{:}, allsubj_mistuning_dev{:}); 
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_montecarlo_not_corected', '.mat')), 'stat');%saved this statistic variable to file

% make the plot
load(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_montecarlo_not_corected', '.mat')));

cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'acticap-64ch-standard2.mat';
cfg.highlight = 'on';
cfg.highlightcolor = 'r';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure;
ft_topoplotER(cfg, grandavg_mistuning_stand)
title('Nonparametric: significant without multiple comparison correction')

%% non-parametric cluster-based
% takes into account that channels correlate with each other, and time as
% well
cfg = [];
cfg.method      = 'triangulation';                         
cfg.layout      = 'acticap-64ch-standard2.mat';                % specify layout of channels
cfg.feedback    = 'yes';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, grandavg_mistuning_stand);

cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [0.2 0.3]; %point by point, still on the predefined time-interval
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 10000;  

Nsub = length(subjects_all);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubj_mistuning_stand{:}, allsubj_mistuning_dev{:});
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_montecarlo_corected_by_neighbours_fixed-window', '.mat')), 'stat');
%stat.posclusters(1); %into command line, to print the results of the biggest cluster
%stat.negclusters(1)

% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'acticap-64ch-standard2.mat';
cfg.highlight = 'on';
cfg.highlightcolor = 'r';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, grandavg_mistuning_stand)
title('Nonparametric: significant with cluster-based multiple comparison correction')

%% cluster-based without predefined window
cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [0 0.7]; %without predefined time-window
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 500;  
cfg.minnbchan        = 2;      % minimal number of neighbouring channels

Nsub = length(subjects_all);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, allsubj_mistuning_stand{:}, allsubj_mistuning_dev{:});
save(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_montecarlo_corected_by_neighbours_nonfixed-window', '.mat')), 'stat');

% make a plot
load(fullfile(path_data, strcat('music_mfeat_', 'mistuning' , '_montecarlo_corected_by_neighbours_nonfixed-window', '.mat')));

cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.']; %1x5 vector, highlight marker symbol series (default ['*', 'x', '+', 'o', '.'] for p < [0.01 0.05 0.1 0.2 0.3]
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
cfg.zlim = [-6 6]; %t-values
cfg.saveaspng = 'cluster';
ft_clusterplot(cfg, stat);

