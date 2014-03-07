function sl06c_realign(cfg, subj)
%SLEEPLIEGE 06C: realign only after fix

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
%-get subject-specific information
%------------%
%-get transform matrix
seg_sn = sprintf('%ss-%s_%04.f_seg_sn.mat', cfg.wmsd, cfg.code, subj);
%------------%

%------------%
%-get functional
rdir = sprintf([cfg.data '%04d/rec/'], subj);
allrdir = dir([rdir 'r*']); % rec folder

for r = 1:numel(allrdir) % r01, r02 etc
  r0dir = [rdir allrdir(r).name filesep];
  allfdir = dir([r0dir 'sf*.nii']);
  
  n_vol = count_volumes_from_name([r0dir allfdir(1).name]);
  for f = 1:n_vol
    sfIMG{r}{f,1} = [r0dir allfdir(1).name ',' num2str(f)];
  end
  
end
%------------%

%---------------------------%
%-do realignment
matlabbatch = [];

matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {seg_sn};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cat(1, sfIMG{:}); % concatenate images
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
  78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 3;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';

spm_jobman('run', matlabbatch)
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

function n_vol = count_volumes_from_name(img_name)

firstvol = str2double(img_name(end-12:end-9));
lastvol = str2double(img_name(end-7:end-4));
n_vol = lastvol - firstvol + 1;