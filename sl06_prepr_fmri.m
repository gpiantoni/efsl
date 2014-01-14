function sl06_prepr_fmri(cfg, subj)
%SLEEPLIEGE 06: PREPR_fMRI

mversion = 7;
%07 11/10/21 feedback on movement parameters
%06 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%05 11/02/11 keep all struct together, in cfg.wmsi, not in rdir
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

%---------------------------%
%-get template images
tpmdir = fullfile(fileparts(which('spm')), 'tpm');
tpm{1,1} = [tpmdir filesep 'grey.nii'];
tpm{2,1} = [tpmdir filesep 'white.nii'];
tpm{3,1} = [tpmdir filesep 'csf.nii'];

fprintf('Running subject %02.f\n', subj);
tic_s = tic;
%---------------------------%

%---------------------------%
%-get structural and functional
sIMG = [];
fIMG = [];
firstIMG = []; % useful to get meanIMG and rp.txt

ndir = sprintf('%04.f', subj);
rdir = [cfg.data, ndir, '/rec/'];

%------------%
%-structural
% stdir = dir([rdir 's*.img']);
% sIMG  = [rdir stdir(1).name];
sIMG   = sprintf('%ss-%s_%04.f.img', cfg.wmsd, cfg.code, subj);
wmsIMG = sprintf('%swms-%s_%04.f.img', cfg.wmsd, cfg.code, subj);

seg_sn = [sIMG(1:end-4) '_seg_sn.mat'];
seg_inv_sn = [sIMG(1:end-4) '_seg_inv_sn.mat'];

%------------%
%-functional
allrdir = dir([rdir 'r*']); % rec folder

for r = 1:numel(allrdir) % r01, r02 etc
  r0dir = [rdir allrdir(r).name filesep];
  allfdir = dir([r0dir 'f*.nii']);
  
  firstIMG(r).dir = r0dir;
  firstIMG(r).img = allfdir(1).name;
  
  n_vol = count_volumes_in_nii([r0dir allfdir(1).name]);
  for f = 1:n_vol
    fIMG{r}{f,1} = [r0dir allfdir(1).name ',' num2str(f)];
  end
  
end
%---------------------------%

%---------------------------%
%-matlabbatch

%------------%
%-Realign
matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estwrite.data = fIMG;

%------%
%-Estimate
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};

%------%
%-Reslice (mean image only)
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch)

%------%
%-Output (mean image)
meanIMG = [firstIMG(1).dir 'mean' firstIMG(1).img ',1'];
%------%

%------%
%-feedback on movement parameters
for r = 1:numel(allrdir) % r01, r02 etc
  rp = [firstIMG(r).dir 'rp_' firstIMG(r).img(1:end-4) '.txt'];
  mp = dlmread(rp); % movement parameters
  totm = sum(abs(diff(mp)));
  mt = sprintf('% 10.4f', totm); % movent text
  output = sprintf('%sr%02.f: %s\n', output, r, mt);
end
%------%
%------------%

%------------%
%-Coreg
matlabbatch = [];

matlabbatch{1}.spm.spatial.coreg.estimate.ref    = {sIMG};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanIMG};
matlabbatch{1}.spm.spatial.coreg.estimate.other  = cat(1, fIMG{:}); % concatenate the cells

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm_jobman('run', matlabbatch)
%------------%

%-----------------%
%-do structural only if necessary
if exist(wmsIMG, 'file')
  
  %-------%
  %-output
  outtmp = sprintf('%s already exists, using it.\n', ...
    wmsIMG);
  output = [output outtmp];
  %-------%
  
else
  
  %-------%
  %-output
  outtmp = sprintf('%s doesn''t exist, it''ll be computed in %s\n', ...
    wmsIMG, cfg.wmsd);
  output = [output outtmp];
  %-------%
  
  %------------%
  %-Segment
  matlabbatch = [];
  
  matlabbatch{1}.spm.spatial.preproc.data = {sIMG};
  
  matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 0];
  matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 0];
  matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
  matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
  matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
  matlabbatch{1}.spm.spatial.preproc.opts.tpm = tpm;
  matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2; 2; 2; 4];
  matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
  matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
  matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
  matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
  matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
  matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
  matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
  
  spm_jobman('run', matlabbatch)
  
  %------%
  %-Output
  % msIMG  = [rdir 'm' stdir(1).name];
  msIMG = sprintf('%sms-%s_%04.f.img', cfg.wmsd, cfg.code, subj);
  
  seg_sn = [sIMG(1:end-4) '_seg_sn.mat'];
  seg_inv_sn = [sIMG(1:end-4) '_seg_inv_sn.mat'];
  %------------%
  
  %------------%
  %-Normalize structural
  matlabbatch = [];
  
  matlabbatch{1}.spm.spatial.normalise.write.subj.matname  = {seg_sn};
  matlabbatch{1}.spm.spatial.normalise.write.subj.resample =  {msIMG};
  matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
  matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
    78 76 85];
  matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
  matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 3;
  matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
  
  spm_jobman('run', matlabbatch)
  
  %------%
  %-Output
  wmsIMG = sprintf('%swms-%s_%04.f.img', cfg.wmsd, cfg.code, subj);
  %-------%
  
  %-------%
  %-cleanup structural
  delete( msIMG)
  delete([msIMG(1:end-3) 'hdr'])
  %-------%
  %------------%
  
end
%-----------------%

if ~cfg.melo
  
  %------------%
  %-Normalize functional
  matlabbatch = [];
  
  matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {seg_sn};
  matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cat(1, fIMG{:}); % concatenate images
  matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
  matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
    78 76 85];
  matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [2 2 2]; % <- only difference with structural
  matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 3;
  matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
  
  spm_jobman('run', matlabbatch)
  
  %------%
  %-Output
  for r = 1:numel(allrdir) % r01, r02 etc
    r0dir = [rdir allrdir(r).name filesep];
    allfdir = dir([r0dir 'wf*.nii']);
    
    n_vol = count_volumes_in_nii([r0dir allfdir(1).name]);
    for f = 1:n_vol
      wfIMG{r}{f,1} = [r0dir allfdir(1).name ',' num2str(f)];
    end
  end
  
  %-------%
  %-cleanup functional
  for r = 1:numel(allrdir)
    delete(fIMG{r}{1}(1:end-2))
    delete([fIMG{r}{1}(1:end-5) 'mat'])
  end
  %------------%

  %------------%
  %-split in pieces because it's too large for some subjects
  for r = 1:numel(wfIMG)
    wfIMG_img{r} = spm_file_split(wfIMG{r}{1}(1:end-2));
    delete(wfIMG{r}{1}(1:end-2))
    
    wfIMG{r} = [];
    for i = 1:numel(wfIMG_img{r})
      wfIMG{r}{i} = [wfIMG_img{r}(i).fname ',1'];
    end
  end

  %------------%
  %-Smooth
  matlabbatch = [];

  matlabbatch{1}.spm.spatial.smooth.data = cat(1, wfIMG{:});

  matlabbatch{1}.spm.spatial.smooth.fwhm = [1 1 1] * cfg.smoo;
  matlabbatch{1}.spm.spatial.smooth.dtype = 0;
  matlabbatch{1}.spm.spatial.smooth.im = 0;
  matlabbatch{1}.spm.spatial.smooth.prefix = 's';

  spm_jobman('run', matlabbatch)
  
  %-------%
  %-cleanup functional
  for r = 1:numel(allrdir)
    for f = 1:numel(wfIMG_img{r})
      delete(wfIMG_img{r}(f).fname)
    end
  end
  
  %------%
  %-Output (merge)
  for r = 1:numel(allrdir) % r01, r02 etc
    r0dir = [rdir allrdir(r).name filesep];
    allfdir = dir([r0dir 'swf*.nii']);
    
    for f = 1:numel(allfdir)
      swfIMG{r}{f,1} = [r0dir allfdir(f).name];
    end
    
    swfIMG_4D = [swfIMG{r}{1}(1:end-10) '.nii'];
    spm_file_merge(swfIMG{r}, swfIMG_4D);
    
    for f = 1:numel(swfIMG{r})
      delete(swfIMG{r}{f})
    end

  end


end

%-------%
%-time
toc_s = toc(tic_s);
fprintf('Analysis finished p%02.f: %s\n', subj, datestr( datenum(0, 0, 0, 0, 0, toc_s), 'HH:MM:SS'));
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

function n_vol = count_volumes_in_nii(img_name)

n_vol = bash(['fslinfo ' img_name ' | awk ''NR==5'' | awk ''{print $2}''']);
n_vol = str2double(n_vol);

function output = bash(command, cwd)

if nargin == 1
  [~, output] = system(['. ~/.bashrc; ' command]);
else
  [~, output] = system(['. ~/.bashrc; cd ' cwd ' ; ' command]);
end