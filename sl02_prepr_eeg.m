function sl02_prepr_eeg(cfg, subj)
%SLEEPLIEGE 02: PREPR_EEG
% althought the batch files are very simple and written by me, it's better
% to use the batch system, even if less transparent.
% If I can change something in the batch, I should get the same results
% from the script and from the batch. Plus, it's easy to create a batch
% from the script if I need it.
% rename_markers.m is our function, no need for the batch there
% 

mversion = 10;
%10 11/09/09 2nd argument for subj (and cfg.subj -> subj) (still obsolete)
%09 11/07/28 obsolete: I ran it three times, I used the best ICA (initially
%            called eeg2). Now that ICA is in recordings/../eeg/conv as the most
%            stable and reliable version. Don't run this again.
%08 11/07/05 use dir struct with raw and preproc dir within each subj dir
%07 11/02/28 workaround: cICA deletes CRC field, we need to put it back
%06 11/02/28 workaround: crc_par.m, l.128, checks for D.other.CRC
%05 11/02/08 uniform cfg
%04 11/02/07 no for subj
%03 11/01/19 bug: we need to set process to 1 every time within loop (and prefix = '')
%02 11/01/18 improved clean-up
%01 11/01/17 created

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

%---------------------------%
%-DONT NOT RUN THIS
output = 'This function is obsolete: Now that ICA is in recordings/../eeg/conv as the most stable and reliable version. Don''t run this again.';
%-----------------%
fprintf(output)
fid = fopen(cfg.log, 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
return
%---------------------------%

%---------------------------%
%-Prepr_eeg (dataset to preprocess)
for i = 1:14
  prepr_eeg{i} = '2';
end

% prepr_eeg{ 2} = '23'; % it's not necessary to preprocess the second (3) dataset
prepr_eeg{14} = '23';
%---------------------------%

for e = 1:numel(prepr_eeg{cfg.subj})
  
  fprintf('Running subject %02.f, dataset %s\n', cfg.subj, prepr_eeg{cfg.subj}(e));
  tic_s = tic;
  
  %---------------------------%
  %-get EEG data
  ndir    = sprintf('%s_%03.f', cfg.code, cfg.subj);
  edir = [cfg.data ndir filesep cfg.eegd filesep];
  alledir = dir([edir '*_' prepr_eeg{cfg.subj}(e) '.vhdr']);
  
  eegdat = alledir(1).name;
  %---------------------------%
  
  %---------------------------%
  %-Read the data
  matlabbatch = [];
  
  matlabbatch{1}.fast.readdata.readbrpr.data = {[edir eegdat]};
  
  spm_jobman('run', matlabbatch)
  
  process = 1;
  prefix = ''; % add prefix as we go
  
  eegdat = [eegdat(1:end-4) 'mat'];
  
  %RENAME
  % eegname = sprintf('%s%s_%03.f_eeg_t1_sleep_r%1.f%s', eeg, cfg.code, cfg.subj, r, ext{e});
  
  %-------%
  %-workaround, bc crc_par.m (l.128) checks whether CRC exists
  load([edir eegdat])
  D.other.CRC = [];
  save([edir eegdat], 'D')
  %---------------------------%
  
  %---------------------------%
  %-Non-linear preprocessing (depending on cfg.ordr)
  while process <= numel(cfg.ordr)
    
    %---------------%
    %-Input data
    
    Dmat = [edir prefix eegdat];
    
    switch cfg.ordr(process)
      
      %------------%
      case 1 % montage
        matlabbatch = [];
        
        matlabbatch{1}.spm.meeg.preproc.montage.D = {Dmat};
        
        matlabbatch{1}.spm.meeg.preproc.montage.montage = {cfg.mont};
        matlabbatch{1}.spm.meeg.preproc.montage.keepothers.yes = 1;
        
        spm_jobman('run', matlabbatch)
        prefix = ['M' prefix];
        
      %------------%
      case 2 % gradient artifact
        
        matlabbatch = [];
        
        matlabbatch{1}.fast.gar.data = {Dmat};
        
        spm_jobman('run', matlabbatch)
        prefix = ['cga_' prefix];
        
      %------------%
      case 3 % pulse artifact
        
        matlabbatch = [];
        
        matlabbatch{1}.fast.pulseartfct.data = {Dmat};
        
        matlabbatch{1}.fast.pulseartfct.options.qrsmethod = 1;
        matlabbatch{1}.fast.pulseartfct.options.bcgmethod = 4;
        matlabbatch{1}.fast.pulseartfct.options.peaksave = 2;
        matlabbatch{1}.fast.pulseartfct.options.nit = 50;
        matlabbatch{1}.fast.pulseartfct.options.ecgchan = 0;
        matlabbatch{1}.fast.pulseartfct.options.badchan = 32;
        
        spm_jobman('run', matlabbatch)
        prefix = ['cICA_' prefix];
        
    end
    
    %------------%
    %-cICA removes "other" field (idk where, complicated code)
    load([edir prefix eegdat])
    if ~isfield(D.other, 'CRC') % <- we created this one just before, if missing, it means it was deleted
      Dold = load(Dmat);
      D.other = Dold.D.other;
      D = meeg(D);
      save(D);
    end
    %------------%
    
    %------------%
    %-delete old data
    delete(Dmat)
    Ddat = [Dmat(1:end-3) 'dat'];
    Deeg = [Dmat(1:end-3) 'eeg'];
    if exist(Ddat, 'file'); delete(Ddat); end
    if exist(Deeg, 'file'); delete(Deeg); end
    %------------%
    
    process = process + 1;
    
  end
  
  %------------%
  %-rename TR markers
  Dmat = [edir prefix eegdat];
  rename_TR_markers(Dmat); % save on disk
  %------------%
  
  %---------------------------%
  %-time
  toc_s = toc(tic_s);
  fprintf('Analysis finished p%02.f, dataset %s: %s\n', cfg.subj, prepr_eeg{cfg.subj}(e), datestr( datenum(0, 0, 0, 0, 0, toc_s), 'HH:MM:SS'));
  %---------------------------%
  
end

%------------%
% delete small files
delete([edir '*.log'])
delete([edir '*.vhdr'])
delete([edir '*.vmrk'])
delete([edir '*.eeg']) % first session
%------------%

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
