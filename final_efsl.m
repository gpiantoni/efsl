function final_efsl(opt)

% 11/12/12 massive restructuring of sl07 and sl09, now comparisons at sl07, ttest or simple f-test at sl09
% 11/12/07 use struct2log instead of custom-made sprintf
% 11/11/29 added wcon to weight contrasts based on number of f2b and b2f slow waves
% 11/10/25 cfg.else instead of cfg.compare
% 11/10/21 use backup markers and skip sl02 and sl03
% 11/10/06 does not depend on savant
% 11/09/28 does not use symbolic link when adding the path
% 11/09/22 include minimum amount of traveling for slow waves
% 11/09/08 ibash was used only across subjects, not use allcond_efsl but loop over condition, running subjects sequentially (slower, but more straightforward)
% 11/09/05 REVERTED: spmA contains date as well (bc the preprocessing might affect analysis)
% 11/09/05 spmA contains date as well (bc the preprocessing might affect analysis)
% 11/09/02 allows second argument, opt, which is merged with cfg (for for-loop in bash)
% 11/09/01 more flexible in specifying contrasts for sl09 (ideally sl08 should use the same scheme)
% 11/08/31 AorB to run main or f2b analysis
% 11/07/28 eegd refers to rec/../conv, which contains eeg2 only
% 11/07/19 wiht ibash, define subjall=1:14 again
% 11/07/15 cfg.evtB is in final_efsl.m (it can be reused)
% 11/07/05 using more consistent file structure
% 11/04/26 always write cfg in log and send email only if not bash for each subject
% 11/04/18 added cfg.dlay
% 11/04/01 log gets unique identifier depending on subjall(1)
% 11/03/31 allow for parallel, across subjects
% 11/03/29 added sl10 and fix cfg.base (actually moved mounting point, not necessary anymore)
% 11/02/11 keep struct in cfg.wmsd folder
% 11/02/09 two analysis: spmA and spmB (slow wave and big/small)
% 11/01/27 all the prepr is here, up to second_level
% 11/01/18 has a log file, send email (needs jvm)
% 11/01/18 created

[~, host] = system('hostname');
if strcmp(host(end-12:end-1), 'partners.org')
  cfg.base = '/PHShome/gp902/projects/efsl/';
  toolboxdir = '/PHShome/gp902/toolbox/';
else
  cfg.base = '/data1/projects/efsl/';
  toolboxdir = '/usr/local/toolbox/';
end

%-------------------------------------%
%-Paths (order is important)----------%
%-------------------------------------%

cfg.scrp = [cfg.base 'scripts/'];

if ~exist('crc_main')
  addpath([cfg.scrp 'final/'])
  addpath([toolboxdir 'FASST_111017/'])
  addpath([toolboxdir 'spm8/'])
  addpath(genpath([toolboxdir 'spm8/external/fieldtrip/']));
  addpath([toolboxdir 'pppi_peak/PPPI/'])
  addpath([matlabroot filesep 'toolbox/stats/stats'])  % nanmean in spm8/fieldtrip is broken
  addpath([toolboxdir 'fieldtrip/qsub/'])
  addpath([toolboxdir 'helpers'])
  
  spm_jobman('initcfg')
  fast = crc_cfg_fasst;
  cfg_util('addapp', fast)
end
spm defaults fmri
%-------------------------------------%

%-------------------------------------%
%-CFG---------------------------------%
%-------------------------------------%

cfg.recd  = [cfg.base 'recordings/efsl/'];
if ~isdir(cfg.recd)
  error('There is no recording directory')
end

cfg.anly = [cfg.base 'analysis/'];
cfg.mdir = [cfg.anly 'masks/'];
cfg.rslt = [cfg.anly 'spm/'];

%-----------------%
%-allow parallel computing, using bash
subjall = [14 8 10 5 11 3 12 7 13 1 9 6 4 2];
cfg.step = [4:13];
HPC = 1;
%-----------------%

%-----------------%
%-cfg sl01_getdata
cfg.code = 'efsl';
cfg.recs = [cfg.recd 'subjects/'];
cfg.data = [cfg.base 'subjects/'];
cfg.rawd = 'raw/';
cfg.type = 'eeg'; % fmri, eeg, both (fmri can be taken directly from rawdata with sl05_divide_rec)
cfg.eegd = 'eeg'; % results can be different depending on ICA (recordings only contains eeg2, the best one)
%-----------------%

%-----------------%
%-cfg sl02_prepr_eeg
% don't run this again, data are already clean
cfg.mont = [cfg.recd 'doc/rerefTP9TP10.mat'];
cfg.ordr = [1 2 3]; % order can be different
%-----------------%

%-----------------%
%-cfg sl03_sw_sp_det
allprefix = {'M', 'cga_', 'cICA_'};
cfg.pref = [allprefix{cfg.ordr(end:-1:1)}]; % you need to add them backwards
cfg.mrkr = [cfg.anly 'fMRI_markers_GP.mat'];
cfg.echk = 'chk';
cfg.hpfilt = 0.2;
cfg.lpfilt = 4;
%-----------------%

%-----------------%
%-cfg sl04_prepare_triggers
%-------%
cfg.swdr = '0'; % '0' or 'range'
cfg.spdr = '0'; % '0' or 'dur'
%-------%

cfg.SWest = 'swstream2'; % 'swstream2'; % or ''
cfg.SWcla = 'mostsouth'; % 'mostsouth' or 'ytype_xylen' or 'beginend' or 'backup' (backup does no have RR)
cfg.mintrvl = 0; % if 'backup', this should be 0
cfg.else = 0; % 1 (f2b) or 2 (b2f) or  0;  if 'backup', this should be ''

cfg.offset = 'no'; % offset due to the rounding error introduced by FASST
%-----------------%

%-----------------%
%-cfg sl05_divide_rec
cfg.minW = 10; % min amount of slow wave to call it a session
cfg.wmsd = [cfg.anly 'avg_mri/']; % don't run seg every time but reuse images
cfg.smoo = 4; % <- bc names change depending on smoothing (although smoothing is applied later)
%-----------------%

%-----------------%
%-cfg sl06_prepr_fmri
cfg.melo = false;
%-----------------%

%-----------------%
%-cfg sl06b_get_melodic
cfg.clme = [cfg.anly 'melodic_feat/clean/'];
cfg.tfsf = [cfg.scrp 'final/private/template.fsf'];  % template fsf
cfg.fixd = [toolboxdir 'fix1.06' filesep];
cfg.ncmp = 20; % number of components
%-----------------%

%-----------------%
%-cfg sl07_first_level
cfg.AorB = 'B'; % 'A' -> only main analysis, 'B' -> only f2b/b2f, 'AB' -> both

cfg.dlay = 0; % delay in seconds
cfg.pmod = ''; % '' or '_pmod-dur' or '_pmod-par' ('_pmod-par' only affects b/sSW contrast, while the main SW contrast remains ''). Note that '_pmod-dur' somehow creates empty columns, so no valid contrasts. idk y, but not important
cfg.bases = 'fir'; % 'hrf' or 'fir'
cfg.basopt = [10 5]; 
% for hrf [0 0] (simple hrf) or [1 1] (with derivative)
% for fir [12 6] length and duration
cfg.volt = 1; % fd1 = no, 2 = yes (big problem if you include it, because there are not enough scans for some subjects!)
cfg.heart = 'yes';
cfg.wcon = 'no'; % yes or no (weight the contrasts for each session based on number of slow waves)

%-------%
%-contrasts to collect (only use one, otherwise it's not correct)
if strcmp(cfg.pmod, '') || strcmp(cfg.pmod, '_pmod-par')
  cfg.evtA(1).img = {'SlowWave'};
  cfg.evtA(1).con  = 1;
  cfg.evtA(2).img = {'Spindle'};
  cfg.evtA(2).con  = 1;
elseif strcmp(cfg.pmod, '_pmod-dur')
  cfg.evtA = {'SW', 'SWpar', 'SP'};
end
%-------%

%-------%
cfg.bsSP = '_bsSP'; % '_bsSP' or '' (include spindles in design matrix of f2b vs b2f)
%-------%

%-------%
%-contrasts to collect
if strcmp(cfg.pmod, '')
  cfg.evtB(1).img = {'f2bSWmain'};
  cfg.evtB(1).con  = 1;
  cfg.evtB(2).img = {'b2fSWmain'};
  cfg.evtB(2).con  = 1;
  cfg.evtB(3).img = {'f2bSWmain' 'b2fSWmain'};
  cfg.evtB(3).con  = [1 -1];
  
elseif strcmp(cfg.pmod, '_pmod-dur') || strcmp(cfg.pmod, '_pmod-par')
  cfg.evtB(1).img = {'f2bSWmain'};
  cfg.evtB(1).con  = 1;
  cfg.evtB(2).img = {'f2bSWparam'};
  cfg.evtB(2).con  = 1;
  cfg.evtB(3).img = {'f2bSWmain' 'f2bSWparam'};
  cfg.evtB(3).con  = [1 1];
  cfg.evtB(3).img = {'f2bSWmain' 'f2bSWparam'};
  cfg.evtB(3).con  = [1 1];
  cfg.evtB(4).img = {'b2fSWmain'};
  cfg.evtB(4).con  = 1;
  cfg.evtB(5).img = {'b2fSWparam'};
  cfg.evtB(5).con  = 1;
  cfg.evtB(6).img = {'b2fSWmain' 'b2fSWparam'};
  cfg.evtB(6).con  = [1 1];
  cfg.evtB(7).img = {'f2bSWparam' 'b2fSWparam'};
  cfg.evtB(7).con  = [1 -1];
end
%-------%
%-----------------%

%-----------------%
%-cfg sl08_second_level_main
cfg.outp = [cfg.anly 'output/'];

%-------%
%-create_mask: mask at the second level (for both sl08 and sl09)
cfg.Dwfu = [toolboxdir 'WFU_PickAtlas_3.0.1/wfu_pickatlas/MNI_atlas_templates/'];
cfg.dMsk = [cfg.anly 'masks/'];

% 'no', 'con1' or other
% I don't trust con1 too much, I think it's only a visual mask, not a p-value mask
% if it's other, then the name you specified is used to create the mask,
% with the ROIs and atlas that you specify. It checks with the current
% cfg.msk2 content matches cfg.atls and ROIs, otherwise, over-writes it.

cfg.msk2 = 'ponsredrawn'; %'ponsredrawn'; %'ponsdrawn'; % 'ponsdrawn' 'lobe_midbrain' 'midbraindrawn' 'ponsdrawn'; % 'no', 'con1' or 'lobe_midbrain' other
% find atlas and regions info:
% /usr/local/toolbox/WFU_PickAtlas_3.0.1/wfu_pickatlas/MNI_atlas_templates
cfg.atls =  'drawn'; %'TD_lobe'; % 'TD_hemisphere' 'TD_lobe' 'atlas71'; % leave empty if you don't use it, or 'drawn'
cfg.ROIs = {''};
%-LC should be in the pons, according to wikipedia
% hemisphere_lfbrainstem is too large
% lobe pons is too anterior, only front part of pons is in the mask
%-------%
%-----------------%

%-----------------%
%-cfg sl09_second_level_SW
%-------%
%-prepare name (don't touch)
for k = 1:numel(cfg.evtA)
  cfg.evtA(k).name = '';
  for i = 1:numel(cfg.evtA(k).img)
    cfg.evtA(k).name = sprintf('%s%s%1.f_', cfg.evtA(k).name, cfg.evtA(k).img{i}, cfg.evtA(k).con(i));
  end
  cfg.evtA(k).name(end) = [];
end

for k = 1:numel(cfg.evtB)
  cfg.evtB(k).name = '';
  for i = 1:numel(cfg.evtB(k).img)
    cfg.evtB(k).name = sprintf('%s%s%1.f_', cfg.evtB(k).name, cfg.evtB(k).img{i}, cfg.evtB(k).con(i));
  end
  cfg.evtB(k).name(end) = [];
end
%-------%
%-----------------%

%-----------------%
%-cfg sl10_roi4ppi
cfg.LCdf = 'act'; % 'ROI_3vxinLC.img'; % 'ROI_...' (if it starts with ROI, it reads it from mask/, with extension) 'act' (activation, FWE, .05) or 'biact' (make act bilateral) maxclu (only largest cluster), bimaxclu (make maxclu bilateral)
cfg.onlyLC = 'no';
cfg.LCic = 1; % index of the contrast of cfg.evtB
cfg.LCpv = 0.05; % p-value to define LC
cfg.LCcr = 'FWE'; % correction to define LC ('FWE' or 'none')
%-----------------%

%-----------------%
%-cfg sl11_ppisubj
cfg.pext = 'eig'; % 'mean' or 'eig'
cfg.pmet = 'trad'; % 'trad' or 'cond' (if 'cond', change ptsk into {'1', 'task_name'}
cfg.ptsk = cfg.LCic ; % which task, index of contrast, similar to cfg.LCic
%-----------------%

%-----------------%
%-cfg sl12_second_level_ppi
%-----------------%

%-----------------%
%-cfg sl13_report
%-------%
%-report results
cfg.gdLC = [6 -38 -22]; % if peak is here, it's LC! (it might depend on voxel size)
cfg.gdLC = [6 -38 -22; -6 -38 -22; 
            6 -38 -24; -6 -38 -24;
            6 -36 -24; -6 -36 -24;
            4 -36 -24; -4 -36 -24]; % if peak is here, it's LC! (it might depend on voxel size)
cfg.csvf = [cfg.anly 'spm/efsl.csv'];
%-------%
%-----------------%
%---------------------------%

%---------------------------%
% OPT -> CFG
%-----------------%
%-merge opt with cfg
if nargin == 1 % opt exists and
  fopt = fieldnames(opt);
  for i = 1:numel(fopt)
    cfg.(fopt{i}) = opt.(fopt{i});
  end
end
%-----------------%

%-----------------%
%-create fields (based on possible changes bc of OPT)
cfg.trigA = sprintf('onsets_sw_dur%s_%03.f-%03.f_sp_dur%s_off%s_%s_%s', ...
  cfg.swdr, cfg.hpfilt*10, cfg.lpfilt*10, cfg.spdr, cfg.offset, cfg.eegd, cfg.pref); % where the trigger is

cfg.trigB = [cfg.trigA '_' cfg.SWest '_' cfg.SWcla sprintf('%03.f', cfg.mintrvl*100) '_' num2str(cfg.else)];

if strcmp(cfg.SWcla, 'backup')
  cfg.trigB = sprintf('sw_dur%s_sp_dur%s_%s', ...
    cfg.swdr, cfg.spdr, cfg.eegd); % where the trigger is
end

if strcmp(cfg.bases, 'hrf');
  firn = ['_hrf' num2str(cfg.basopt(1)), num2str(cfg.basopt(2))]; % hrf and two derivatives
else
  firn = sprintf('_l%02.fo%02.f', cfg.basopt(1), cfg.basopt(2));
end
cfg.spmA = ['spm_s' num2str(cfg.smoo) '_' cfg.pref '_sw' cfg.swdr cfg.pmod sprintf('_d%02.f', cfg.dlay*10) firn '_sp' cfg.spdr '_' cfg.eegd '_v' num2str(cfg.volt) '_h' cfg.heart '_w' cfg.wcon];
cfg.dirA = [cfg.rslt cfg.spmA filesep];

cfg.spmB = [cfg.spmA '_' cfg.SWest '_' cfg.SWcla sprintf('%03.f', cfg.mintrvl*100)  '_' num2str(cfg.else) cfg.bsSP];
cfg.dirB = [cfg.rslt cfg.spmB filesep];
%-----------------%
%---------------------------%

%---------------------------%
%-Log file
cfg.log = sprintf('%slog/prepr_%s_%s_%1.f.txt', ...
  cfg.anly, datestr(now, 'yy-mm-dd'), datestr(now, 'HH_MM_SS'), subjall(1));
fid = fopen(cfg.log, 'w');
output = sprintf('Analysis started at %s on %s\n', ...
  datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
fprintf(output)
fwrite(fid, output);
[~, gitversion] = system(['git rev-parse HEAD']);
fwrite(fid, sprintf('GIT VERSION: %s\n', gitversion));

%-----------------%
%-cfg in log
output = struct2log(cfg);
fprintf(output)
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-SINGLE-SUBJECT ANALYSIS-------------%
%-------------------------------------%
subjcell = num2cell(subjall);
cfgcell = repmat({cfg}, 1, numel(subjall));
rmdir([cfg.scrp 'final/qsublog'], 's')
mkdir([cfg.scrp 'final/qsublog'])
cd([cfg.scrp 'final/qsublog'])

%---------------------------%
%-if running the exact same analysis, remove previous results
if any(cfg.step ==  7) && numel(subjall) == 14
  if isdir(cfg.dirA); rmdir(cfg.dirA, 's'); mkdir(cfg.dirA); end
  if isdir(cfg.dirB); rmdir(cfg.dirB, 's'); mkdir(cfg.dirB); end
end
%---------------------------%

%---------------------------%
%-Copy original data
% it deletes subject dir!
if any(cfg.step ==  1)
  disp('running sl01_getdata')
  if HPC
    qsubcellfun(@sl01_getdata, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
  else
    sl01_getdata(cfgcell{1}, subjcell{1})
  end
end
%---------------------------%

%---------------------------%
%-EEG preprocessing
if any(cfg.step ==  2) && 0
  
  if strcmp(cfg.SWcla, 'backup')
    warning('using backup markers: sl02_prepr_eeg is not necessary')
  else
    disp('running sl02_prepr_eeg')
    qsubcellfun(@sl02_prepr_eeg, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
  end
  
end
%---------------------------%

%---------------------------%
%-slow wave and spindle detection
if any(cfg.step ==  3)
  if strcmp(cfg.SWcla, 'backup')
    warning('using backup markers: sl03_sw_sp_det is not necessary')
  else
    disp('running sl03_sw_sp_det')
    if HPC
      qsubcellfun(@sl03_sw_sp_det, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
    else
      sl03_sw_sp_det(cfgcell{1}, subjcell{1})
    end
  end
end
%---------------------------%

%---------------------------%
%-prepare triggers
if any(cfg.step ==  4)
  disp('running sl04_prepare_triggers')
  if HPC
    qsubcellfun(@sl04_prepare_triggers, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
  else
    sl04_prepare_triggers(cfgcell{1}, subjcell{1})
  end
end
%---------------------------%

%---------------------------%
%-Create recording folders
% this should depend on prepare_triggers as well, one period w/out SW shouldn't be used
if any(cfg.step ==  5)
  disp('running sl05_divide_rec')
  if HPC && false  % sshfs is only mounted on the pc running matlab
    qsubcellfun(@sl05_divide_rec, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
  else
    for i = 1:numel(subjcell)
      sl05_divide_rec(cfgcell{i}, subjcell{i})
    end
  end
end
%---------------------------%

%---------------------------%
%-fMRI preprocessing
if any(cfg.step ==  6)
   disp('running sl06_prepr_fmri')
    if HPC 
      qsubcellfun(@sl06_prepr_fmri, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
    else
      sl06_prepr_fmri(cfgcell{1}, subjcell{1})
    end
  if cfg.melo
     disp('running sl06b_get_melodic')
    if HPC
      qsubcellfun(@sl06b_run_melodic, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab');
    else
      sl06b_run_melodic(cfgcell{1}, subjcell{1})
    end
  end    
end
%---------------------------%

%---------------------------%
%-first-level stats
% don't remove tdir for the moment, include parametric!
if any(cfg.step ==  7)
  disp('running sl07_first_level')
  if HPC
    qsubcellfun(@sl07_first_level, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab');
  else
    sl07_first_level(cfgcell{1}, subjcell{1})
  end
end
%---------------------------%
%-------------------------------------%

%-------------------------------------%
%-SECOND-LEVEL RESULTS----------------%
%-------------------------------------%

%---------------------------%
%-second-level main
if any(cfg.step ==  8) && ~isempty(strfind(cfg.AorB, 'A'))
  create_masks(cfg)
  disp('running sl08_second_level_main')
  sl08_second_level_main(cfg)
end
%---------------------------%

%---------------------------%
%-second-level big vs small
if any(cfg.step ==  9) && ~isempty(strfind(cfg.AorB, 'B'))
  create_masks(cfg)
  disp('running sl09_second_level_SW')
  sl09_second_level_SW(cfg)
end
%---------------------------%

%---------------------------%
%-create roi for PPI (LC)
if any(cfg.step == 10) && ~strcmp(cfg.LCdf(1:3), 'ROI')
  disp('running sl10_roi4ppi')
  sl10_roi4ppi(cfg)
end
%---------------------------%

if strcmp(cfg.LCdf(1:3), 'ROI') || exist([cfg.dirB cfg.evtB(cfg.LCic).name filesep 'ROI_' cfg.LCdf '.hdr'], 'file')
  
  %---------------------------%
  %-if running the exact same analysis, remove previous results
  if any(cfg.step ==  11) && numel(subjall) == 14
    dirP = [cfg.dirB cfg.evtB(cfg.LCic).name filesep 'PPI_' cfg.LCdf filesep];
    if isdir(dirP); rmdir(dirP, 's'); mkdir(dirP); end
  end
  %---------------------------%
  
  %---------------------------%
  %-ppi for each subject
  if any(cfg.step == 11)
    cd([cfg.scrp 'final/qsublog'])
    disp('running sl11_ppi_subj')
    if HPC
      qsubcellfun(@sl11_ppisubj, cfgcell, subjcell, 'memreq', [], 'timreq', [], 'queue', 'matlab')
    else
      sl11_ppisubj(cfgcell{1}, subjcell{1})
    end
  end
  %---------------------------%
  
  %---------------------------%
  %-ppi across subjects
  if any(cfg.step == 12)
    disp('running sl12_second_level_ppi')
    sl12_second_level_ppi(cfg)
  end
  %---------------------------%
  
end

%---------------------------%
%-ppi across subjects
if any(cfg.step == 13)
  sl13_report(cfg)
end
%---------------------------%

%-----------------%
%-send email
% send_email(cfg)
%-----------------%

cd([cfg.scrp 'final/'])
%-------------------------------------%
