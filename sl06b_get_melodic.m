function sl06b_get_melodic(cfg, subj)
%SLEEPLIEGE 06B: get fMRI data cleaned up with Melodic
% takes output of melodic in the analysis folder and recreates 'swf' files,
% with extension .hdr/.img using the correct name

mversion = 1;

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
%-loop over recordings
rdir = sprintf([cfg.data '%04d/rec/'], subj);
rsess = dir([rdir 'r*']);

for i_r = 1:numel(rsess)
  
  %-----------------%
  recdir = [rdir rsess(i_r).name filesep];
  meldir = [recdir 'melodic' filesep];
  if isdir(meldir)
    rmdir(meldir, 's')
  end
  mkdir(meldir)
  %-----------------%
  
  clean_niigz = sprintf('%sf%02d-%s_s%d_nl.nii.gz', cfg.clme, subj, rsess(i_r).name, cfg.smoo);
  mel_niigz = sprintf('%sf%02d-%s_s%d.nii.gz', meldir, subj, rsess(i_r).name, cfg.smoo);
  copyfile(clean_niigz, mel_niigz);
  bash(['fslsplit ' mel_niigz ' wf'], meldir);
  delete(mel_niigz)
  
  % old files, based on order  
  f_mri = dir([recdir 'f*.hdr']);  % use f at the beginning to avoid swf files already present
  f_mri_sorted = sort({f_mri.name});  % to be safe
  
  % new files, whose name should match the names of the old files
  m_mri = dir([meldir '*.nii.gz']);
  m_mri_sorted = sort({m_mri.name});  % to be safe
  
  for i_m = 1:numel(f_mri_sorted) 
    niigz = [meldir m_mri_sorted{i_m}];
    hdrimg = [recdir 'sw' f_mri_sorted{i_m}];
    bash(['fslchfiletype NIFTI_PAIR ' niigz ' ' hdrimg]);
  end
  
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

function output = bash(command, cwd)

if nargin == 1
  [~, output] = system(['. ~/.bashrc; ' command]);
else
  [~, output] = system(['. ~/.bashrc; cd ' cwd ' ; ' command]);
end


