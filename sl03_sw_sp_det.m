function sl03_sw_sp_det(cfg, subj)
%SLEEPLIEGE 03: Slow wave, spindle detection

mversion = 10;
%10 14/03/05 use subset of functions of fast
%09 11/09/28 change filter properties to detect sw
%08 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%07 11/07/28 get the right name
%06 11/07/05 edir is within efsl_00XX/preproc/eeg
%05 11/02/23 added scorer field for crc_cfg_sw
%04 11/02/08 uniform cfg
%03 11/02/07 no for subj
%02 11/01/16 using matlab batch (instead of scripting)
%01 11/01/14 created

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

%------------%
% divide fMRI scans into sessions
load(cfg.mrkr)

fprintf('Running subject %02.f\n', subj);
tic_s = tic;
%------------%

%---------------------------%
%-mkr is complicated for p02 and p14
if subj == 2
  if any(mkr(2).mkr(:) > 535) % all the mkr should be below 535, end of first file
    error('check markers for p02, probl it needs to be implemented')
  else
    mkrS = mkr(subj).mkr; % normal
  end
  
elseif  subj == 14
  outbound_mkr = mkr(14).mkr(:,1) > 2629;
  
  mkrS = mkr(subj).mkr( ~outbound_mkr, :);
  
else
  mkrS = mkr(subj).mkr;
end
%---------------------------%

%---------------------------%
%-get EEG data
edir = sprintf('%s%04.f/%s/', cfg.data, subj, cfg.eegd);

eegdat = sprintf('%s_%04.f_%s_sleep_2.mat', ...
  cfg.code, subj, cfg.eegd);
%---------------------------%

%---------------------------%
%-make chunks
matlabbatch = [];

matlabbatch{1}.fasst.chunking.data = {[edir eegdat]};

for c = 1:size(mkrS,1)
  matlabbatch{1}.fasst.chunking.chunk(c).chunk_beg.t_mark.m_type = mkrS(c,1);
  matlabbatch{1}.fasst.chunking.chunk(c).chunk_beg.t_mark.m_ind = 1;
  matlabbatch{1}.fasst.chunking.chunk(c).chunk_end.t_mark.m_type = mkrS(c,2);
  matlabbatch{1}.fasst.chunking.chunk(c).chunk_end.t_mark.m_ind = 1;
end

matlabbatch{1}.fasst.chunking.options.overwr = 1;
matlabbatch{1}.fasst.chunking.options.fn_prefix = cfg.echk;
matlabbatch{1}.fasst.chunking.options.numchunk = 1;

spm_jobman('run', matlabbatch)

%-----------------%
%-subject 14 has two sessions in EEG (this is the second)
if  subj == 14 && any(outbound_mkr)
  
  eegdat3 = sprintf('%s_%04.f_%s_sleep_3.mat', ...
    cfg.code, subj, cfg.eegd);

  mkrS = mkr(subj).mkr( outbound_mkr, :);
  
  %-------%
  matlabbatch = [];
  
  matlabbatch{1}.fasst.chunking.data = {[edir eegdat3]};
  
  for c = 1:size(mkrS,1)
    matlabbatch{1}.fasst.chunking.chunk(c).chunk_beg.t_mark.m_type = mkrS(c,1);
    matlabbatch{1}.fasst.chunking.chunk(c).chunk_beg.t_mark.m_ind = 1;
    matlabbatch{1}.fasst.chunking.chunk(c).chunk_end.t_mark.m_type = mkrS(c,2);
    matlabbatch{1}.fasst.chunking.chunk(c).chunk_end.t_mark.m_ind = 1;
  end
  
  matlabbatch{1}.fasst.chunking.options.overwr = 1;
  matlabbatch{1}.fasst.chunking.options.fn_prefix = cfg.echk;
  matlabbatch{1}.fasst.chunking.options.numchunk = numel(find(outbound_mkr == 0)) + 1; % start counting from the first chunk AFTER the chunks from the first EEG file
  
  spm_jobman('run', matlabbatch)
  %-------%
end
%-----------------%
%---------------------------%

%---------------------------%
%-loops through chunks
for c = 1:size(mkr(subj).mkr,1)
  
  %-----------------%
  %-get chunk data
  if  subj == 14 && mkr(subj).mkr(c,1) > 2629
    chkdata = [edir cfg.echk num2str(c) '_' eegdat3];
  else
    chkdata = [edir cfg.echk num2str(c) '_' eegdat];
  end
  %-----------------%
  
  %-----------------%
  %-detect slow waves
  handles = [];
  handles.fname = chkdata;
  handles.highfc = cfg.hpfilt;
  handles.lowfc = cfg.lpfilt;
 
  handles.analyse = 2;  % whole file
  handles.fmri = false;
  handles.reref = true;
  handles.roisel = true;
  handles.review = false;
  
  crc_SWS_detect(handles)
  %-----------------%
  
  %-----------------%
  %-detect spindles
  matlabbatch = [];
  
  matlabbatch{1}.fasst.wavedetect.sp.data = {chkdata};
  
  matlabbatch{1}.fasst.wavedetect.sp.sel.allf = true;
  matlabbatch{1}.fasst.wavedetect.sp.reref = true;
  matlabbatch{1}.fasst.wavedetect.sp.filt.hpfilt = 11;
  matlabbatch{1}.fasst.wavedetect.sp.filt.lpfilt = 20;
  matlabbatch{1}.fasst.wavedetect.sp.review = false;
  matlabbatch{1}.fasst.wavedetect.sp.wavlet = false;
  spm_jobman('run', matlabbatch)
  %-----------------%
  
end
%---------------------------%

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
