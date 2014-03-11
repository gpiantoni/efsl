function sl06b_run_featfix(cfg, subj)
%SLEEPLIEGE 06B: run featfix to clean up fMRI

mversion = 3;

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
  %-create directory
  recdir = [rdir rsess(i_r).name filesep];
  fixdir = [recdir 'fix.feat' filesep];
  if isdir(fixdir)
    rmdir(fixdir, 's')
  end
  
  design_fsf = [recdir 'design.fsf'];
  %-----------------%
  
  %-----------------%
  %-create template
  f = fopen(cfg.tfsf, 'r');
  template = fscanf(f, '%c'); 
  fclose(f);

  fsf = strrep(template, 'XXX_outputdir', fixdir(1:end-6));
  fsf = strrep(fsf, 'XXX_smooth', num2str(cfg.smoo));

  pattern = dir([recdir sprintf('f%02d-r%02d-s*-*-*.nii', subj, i_r)]);
  dirty_fmri = [recdir pattern(1).name];

  [dirname, filename, ext] = fileparts(dirty_fmri);
  fsf = strrep(fsf, 'XXX_feat_files', fullfile(dirname, filename));
  npts = bash(['fslinfo ' dirty_fmri ' | gawk ''FNR == 5 {print $2}'' ']);
  fsf = strrep(fsf, 'XXX_npts', strtrim(npts));
  sIMG = sprintf('%ss-%s_%04d', cfg.wmsd, cfg.code, subj);  % no ext
  sIMG_brain = sprintf('%ss-%s_%04d_brain', cfg.wmsd, cfg.code, subj);  % bet'ed
  if ~exist(sIMG_brain)
    bash(['bet ' sIMG ' ' sIMG_brain]);
  end
  fsf = strrep(fsf, 'XXX_highres_files', sIMG_brain);
  
  f = fopen(design_fsf, 'w');
  fprintf(f, strrep(fsf, '%', '%%'));
  fclose(f);
  %-----------------%

  %-----------------%
  %-run feat
  output = [output sprintf('feat started at %s\n', datestr(now, 'HH:MM:SS'))];
  bash(['feat ' design_fsf]);  
  output = [output sprintf('feat ended at %s\n', datestr(now, 'HH:MM:SS'))];
  %-----------------%
  
  %-----------------%
  %-run fix
  output = [output sprintf('fix started at %s\n', datestr(now, 'HH:MM:SS'))];
  bash([cfg.fixd 'fix ' fixdir(1:end-1) ' ' cfg.fixd 'training_files/Standard.RData ' num2str(cfg.ncmp)])
  output = [output sprintf('fix ended at %s\n', datestr(now, 'HH:MM:SS'))];
  %-----------------%

  %-----------------%
  %-move file as unzipped to main folder
  clean_fmri = [fixdir 'filtered_func_data_clean.nii.gz'];
  ready_fmri = fullfile(dirname, ['s' filename, '.nii.gz']);
  copyfile(clean_fmri, ready_fmri)
  
  gunzip(ready_fmri)
  delete(ready_fmri)
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

function output = bash(command, cwd)

% Use system libraries
MatlabPath = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', getenv('PATH'))

if nargin == 1
  [~, output] = system(['. ~/.bashrc; ' command]);
else
  [~, output] = system(['. ~/.bashrc; cd ' cwd ' ; ' command]);
end

setenv('LD_LIBRARY_PATH', MatlabPath)
