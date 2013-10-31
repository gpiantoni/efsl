function sl04_prepare_triggers(cfg, subj)
%SLEEPLIEGE 04: PREPARE_TRIGGERS

mversion = 15;
%15 11/12/08 there is a small difference between the eeg marker and fmri volume
%14 11/12/07 include RR detection
%13 11/10/21 don't compute them if using backup
%12 11/09/22 include minimum amount of traveling for slow waves
%11 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%10 11/07/05 edir is within efsl_00XX/eeg
%09 11/03/07 add SW_param for big and small waves
%08 11/02/10 extra output
%07 11/02/09 no waves: 0, big wave: 1, small wave: 2
%06 11/02/08 uniform cfg
%05 11/02/08 added classify_SW
%04 11/02/07 no for subj
%03 11/01/27 allows for SW_dur = range
%02 11/01/17 it saves SP_dur (but it's zeros, like SW_dur)
%01 11/01/17 created

% but if there are no slow waves, you have to remove the sessions
% only from fMRI model specification, delete the sessions without SW

%-----------------%
%-input
if nargin == 1
  subj = cfg.subj;
end
%-----------------%

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s (v%02.f) started at %s on %s\n', ...
  subj, mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-prepare directories
ndir = sprintf('%04.f', subj);
trdir = sprintf('%s%s%s', cfg.data, ndir, '/spm/triggers/');
if ~exist(trdir, 'dir')
  mkdir(trdir)
end
%-----------------%
%-where the backup directory is
backupdir = '/data/projects/efsl/backup/projects/efsl/subjects/';
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-use BACKUP
if strcmp(cfg.SWcla, 'backup')
  
  %---------------------------%
  %-main SW
  %--------%
  trdirbkp  = sprintf('%s%04.f%s', backupdir, subj, '/spm/triggers/');
  filebkp = [trdirbkp 'onsets_sw_dur0_sp_dur0_eeg_cICA_cga_M.mat'];
  %--------%
  
  copyfile(filebkp, [trdir cfg.trigA '.mat'])
  %---------------------------%
  
  %---------------------------%
  %-big vs small
  %--------%
  trdirbkp  = sprintf('%s%04.f%s', backupdir, subj, '/spm/triggers/');
  filebkp = [trdirbkp 'onsets_sw_dur0_sp_dur0_eeg_cICA_cga_M_swstream2_mostsouth.mat'];
  %--------%
  
  copyfile(filebkp, [trdir cfg.trigB '.mat'])
  %---------------------------%
  
  return
  
end
%-------------------------------------%

%-------------------------------------%
%-directory for main analysis
if exist([trdir 'temp'], 'dir')
  rmdir([trdir 'temp'], 's')
end

edir = [cfg.data ndir filesep cfg.eegd filesep];
pattern = [edir 'chk*.mat'];
chks = dir(pattern);

if isempty(chks)
  disp(['no chunks matching ' cfg.pref ' in ' edir])
  return
end
%-------------------------------------%

%-------------------------------------%
%-main effect triggers
SW_onset = [];
SW_dur   = [];
SP_onset = [];
SP_dur   = [];
RR = [];
load(cfg.mrkr) % used by R-R

for r = 1:numel(chks)
  load([edir chks(r).name])
  
  %-----------------%
  %-measure R-R
  ECGchan = find(strcmp({D.channels.label}, 'ECG'));
  
  peaks = fmrib_qrsdetect(meeg(D), ECGchan, 0); % Detect peaks, in samples over recording
  
  heart  = zeros(1, D.Nsamples);
  heart(peaks) = [diff(peaks/D.Fsample) 1]; % if an R was present in a particular sample

  %-------%
  %-find number of volume
  %TODO: why is the n of volumes different between eeg and fmri (problem in chunking tool)
  RT = 2.46; % SPM.xY.RT
  
  RTs = RT * D.Fsample;
  nvoleeg = ceil(D.Nsamples / RTs);
  
  nvolfmri = diff(mkr(subj).mkr(r,:))+1;
  %-------%
  
  beats = NaN(nvolfmri, 1);
  
  for i = 1:nvoleeg
    
    begvol = RTs * (i-1) + 1;
    endvol = RTs * i;
    if endvol > D.Nsamples
      endvol = D.Nsamples;
    end
    
    hvol = heart(begvol:endvol);
    beats(i) = mean(hvol(find(hvol)));
    
  end
  
  beats(isnan(beats)) = nanmean(beats);
  mbeats = mean(beats);
  beats = beats - mbeats;
  
  RR{r} = beats(1:nvolfmri); % in some cases, nvoleeg > nvolfmri !
  
  %-------%
  %-output
  outtmp = sprintf('r%1.f: number of volumes for EEG is % 4.f and for fMRI is % 4.f\n   ECG index is% 3.f, mean heart rate is %1.3f\n', ...
    r, nvoleeg, nvolfmri, ECGchan, mbeats);
  output = [output outtmp];
  %-------%
  %-----------------%
  
  %-----------------%
  %-offset due to marker shift in FASST
  begmkr = find([D.trials.events.value] == mkr(subj).mkr(r,1));
  if strcmp(cfg.offset, 'yes')
    offset = D.trials.events(begmkr).time;
  else
    offset = 0;
  end
  
  %-------%
  %-output
  outtmp = sprintf('r%1.f: marker% 4.f is event n. % 2.f with offset %5.3fs\n', ...
    r, mkr(subj).mkr(r,1), begmkr, offset);
  output = [output outtmp];
  %-------%
  %-----------------%
  
  D = meeg(D);
  if isfield(D.CRC, 'SW')
    for w = 1:numel(D.CRC.SW.SW)
      SW_onset{r}(w,1) = D.CRC.SW.SW(w).negmax_tp / fsample(D) - offset;
      
      if strcmp(cfg.swdr, '0')
        SW_dur{r}(w,1)   = 0;
      elseif strcmp(cfg.swdr, 'range')
        SW_dur{r}(w,1)   = range(D.CRC.SW.SW(w).delays)/ 1000; % convert ms -> s
      end
      
    end
  else
    disp(['for subj ' num2str(subj) ' chunk ' num2str(r) ' doesn''t have SW'])
    SW_onset{r} = [];
    SW_dur{r} = [];
  end
  
  if isfield(D.CRC, 'spindles')
    SP_onset{r}(:,1)= D.CRC.spindles.bounds(:,1) / fsample(D) - offset;
    if strcmp(cfg.spdr, '0');
      SP_dur{r}       = zeros(size(SP_onset{r}));
    elseif strcmp(cfg.spdr, 'dur');
      SP_dur{r}       = D.CRC.spindles.duration / 1000; % convert ms -> s
    end
  else
    disp(['for subj ' num2str(subj) ' chunk ' num2str(r) ' doesn''t have spindles'])
    SP_onset{r} = [];
    SP_dur{r} = [];
  end
  
  %-------%
  %-output
  outtmp = sprintf('r%1.f: number of SW %1.f, number of SP %1.f\n', ...
    r, numel(SW_onset{r}), numel(SP_onset{r}));
  output = [output outtmp];
  %-------%
  
end

%-----------------%
%-save main effect triggers
disp(['subj ' num2str(subj) ' has correct markers'])
save([trdir cfg.trigA], 'SW_onset', 'SW_dur', 'SP_onset', 'SP_dur', 'RR')
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-divide big vs small SW
if ~isempty(cfg.SWest)
  
  %---------------------------%
  %-recalculate them
  bSW_onset = [];
  sSW_onset = [];
  bSW_dur = [];
  sSW_dur = [];
  
  for r = 1:numel(chks)
    
    %-----------------%
    %-classify SW
    load([edir chks(r).name])
    D = meeg(D);
    [SWtype, SWparam] = classify_SW(cfg, D);
    %-----------------%
    
    bSW_onset{r} = SW_onset{r}(SWtype == 1);
    bSW_dur{r}   = SW_dur{r}(SWtype == 1);
    bSW_param{r} = SWparam(SWtype == 1);
    
    sSW_onset{r} = SW_onset{r}(SWtype == 2);
    sSW_dur{r}   = SW_dur{r}(SWtype == 2);
    sSW_param{r} = SWparam(SWtype == 2);
    
    %-------%
    %-output
    outtmp = sprintf('r%1.f: number of bSW %1.f, number of sSW %1.f\n', ...
      r, numel(bSW_onset{r}), numel(sSW_onset{r}));
    output = [output outtmp];
    %-------%
    
  end
  
  %-----------------%
  %-save main effect triggers
  save([trdir cfg.trigB], 'bSW_onset', 'bSW_dur', 'bSW_param', 'sSW_onset', 'sSW_dur', 'sSW_param')
  clear bSW* sSW*
  %-----------------%
  %---------------------------%
end
%-------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s (v%02.f) ended at %s on %s after %s\n\n', ...
  subj, mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen(cfg.log, 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
