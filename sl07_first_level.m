function sl07_first_level(cfg, subj)
%SLEEPLIEGE 07: FIRST LEVEL

mversion = 13;
%13 11/12/12 massive reconstruction, constrasts are computed here, second-level is only for difference from the mean
%12 11/12/07 include RR as regressor of no-interest
%11 11/11/29 added wcon to weight contrasts based on number of f2b and b2f slow waves
%10 11/09/22 include minimum amount of traveling for slow waves
%09 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%08 11/04/22 more consistent names for contrasts
%07 11/04/18 added cfg.dlay
%06 11/04/13 either main or contrast of interest (faster)
%05 11/02/16 allows for FIR, flexible contrasts. Double-check that numbering of contrasts matches the conditions
%04 11/02/08 uniform cfg
%03 11/02/07 no for subj
%02 11/01/18 fixed small bugs, copy contrasts
%01 11/01/17 created

spm_jobman('initcfg');

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
%-get directories, triggers, img
%-------%
%-get directories
ndir = sprintf('%04.f', subj);
rdir = sprintf('%s%s/rec/', cfg.data, ndir);
sdir = sprintf('%s%s/spm/', cfg.data, ndir);
%-------%

%------%
%-swf img

%-use swf (spm-only pipeline) or wsf (feat/fix pipeline)
if cfg.melo
  preproc = 'wsf';
else
  preproc = 'swf';
end  

allrdir = dir([rdir 'r*']); % rec folder

clear swfIMG rp

for r = 1:numel(allrdir) % r01, r02 etc
  r0dir = [rdir allrdir(r).name filesep];
  allfdir = dir([r0dir preproc '*.nii']);
  
  n_vol = count_volumes_from_name([r0dir allfdir(1).name]);
  for f = 1:n_vol
    swfIMG{r}{f,1} = [r0dir allfdir(1).name ',' num2str(f)];
  end
  
  swfIMG_3D{r} = split_into_3D(swfIMG{r}{1});
  swfIMG{r} = swfIMG_3D{r};
  
  % D has also regressed out any residual effects of movement using the
  % movement parameters from affine realignment of the volumes, don't do
  % any further realignment in SPM (and of course don't include motion
  % parameters in your stats model in SPM)
  if cfg.cvrp
    allrp = dir([r0dir 'rp*.txt']);
    rp{r} = [r0dir allrp(1).name];
  end
  
end
%-------%

%-------%
%-triggers
trigfile = [sdir 'triggers/' cfg.trigA '.mat'];
if exist(trigfile, 'file')
  load(trigfile)
else
  error('could not find trigger file')
end
%-------%

%-------%
%-triggers
% these are the triggers for the big/small comparison. However, i need to
% load them here, instead of after the if ~isempty(cfg.SWest) because these
% triggers are used to check whether a session should be run or not (even
% in the main analysis). It doesn't make sense to run the main analysis on
% one dataset and the big/small contrasts on slighly different data.

trigfile = [sdir 'triggers/' cfg.trigB '.mat'];
if exist(trigfile, 'file')
  load(trigfile)
else
  error('could not find f2b/b2f SW trigger file')
end
%-------%
%-------------------------------------%

%-------------------------------------%
%-Specify first-level model: MAIN
if strfind(cfg.AorB, 'A')
  
  %-------%
  %-glm folder
  gdir = [sdir 'glm' filesep];
  
  if exist(gdir, 'dir'); rmdir(gdir, 's'); end
  mkdir(gdir)
  %-------%
  
  %-------%
  %-contrast folder
  conA = [cfg.dirA 'contrasts' filesep];
  
  if ~exist(conA, 'dir'); mkdir(conA); end
  ext = {'hdr', 'img'};
  %-------%
  
  %-----------------%
  %-specify model
  matlabbatch = [];
  
  %-------%
  %-first level parameters
  matlabbatch{1}.spm.stats.fmri_spec.dir = {gdir};
  matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
  matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.46;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
  %-------%
  
  %-------%
  %-loop through recordings (sessions)
  % use two counters: r and rcnt
  % r is the number of sessions, indep of n of SWs (use that for SW_onset and SW_dur)
  % rcnt depends on cfg.minW (does not include sessions with too few slow waves). Use that one for sess, swfIMG and rp
  rcnt = 0;
  rall = [];
  for r = 1:numel(SW_onset)
    
    if size(SW_onset{r},1) >= cfg.minW && ...
        size(bSW_onset{r},1) >= cfg.minW && ...
        size(sSW_onset{r},1) >= cfg.minW
      rcnt = rcnt + 1; % <- rcnt is the sessions actually written on file (enough SW), r is the sessions from sleep scoring
      rall = [rall r]; % <- we used it when writing contrasts to know which sessions are actually used
      
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).scans = swfIMG{rcnt};
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).name = 'SlowWave';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).onset = SW_onset{r}' + cfg.dlay;
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).duration = SW_dur{r}';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).tmod = 0;
      if strcmp(cfg.pmod, '') || strcmp(cfg.pmod, '_pmod-par')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
      elseif strcmp(cfg.pmod, '_pmod-dur')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).name  = 'dur';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).param = SW_dur{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).poly  = 1;
      end
      
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).name = 'Spindle';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).onset = SP_onset{r};
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).duration = SP_dur{r};
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).tmod = 0;
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
      
      if strcmp(cfg.heart, 'yes')
         matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress.name = 'heart';
         matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress.val  = RR{r};
      else
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress = struct('name', {}, 'val', {});
      end
      
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).multi = {''};
      if cfg.cvrp
         matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).multi_reg = rp(rcnt);
      end
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).hpf = 128;
    end
    
  end
  %-------%
  
  %-------%
  %-bases
  if strcmp(cfg.bases, 'hrf')
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = cfg.basopt;
  else
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = cfg.basopt(1);
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order  = cfg.basopt(2);
  end
  %-------%
  
  %-------%
  %-first level parameters
  matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
  matlabbatch{1}.spm.stats.fmri_spec.volt = cfg.volt;
  matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
  matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  %-------%
  
  spm_jobman('run', matlabbatch)
  
  %---------------------------%
  %-Estimate first-level model
  matlabbatch = [];
  
  matlabbatch{1}.spm.stats.fmri_est.spmmat = {[gdir 'SPM.mat']};
  matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
  
  spm_jobman('run', matlabbatch)
  %---------------------------%
  
  %---------------------------%
  %-Create contrasts  

  load([gdir 'SPM.mat'])
  
  %-----------------%
  %-prepare contrasts
  %-------%
  %-contrast weight (for each session)
  if strcmp(cfg.wcon, 'yes')
    wcon = []; % weight for the contrast
    wcon{1} = cellfun(@numel, SW_onset);
    wcon{1} = wcon{1}(rall); % it's important to only keep the sessions which are used (have enough events)
    wcon{1} = wcon{1} / sum(wcon{1}); % rescale to one
    
    wcon{2} = cellfun(@numel, SP_onset);
    wcon{2} = wcon{2}(rall); % it's important to only keep the sessions which are used (have enough events)
    wcon{2} = wcon{2} / sum(wcon{2}); % rescale to one
    
  else
    wcon{1} = ones(1, numel(rall)) / numel(rall);
    wcon{2} = ones(1, numel(rall)) / numel(rall);
  end
  wcon{3} = ones(1, numel(rall)) / numel(rall); % for all the other cases
  %-------%
  
  %-------%
  %-number of bases, contrasts, betas
  if strcmp(cfg.bases, 'hrf');
    nbases = sum(cfg.basopt) + 1;
  else
    nbases =  cfg.basopt(2);
  end
  ncon = numel(cfg.evtA);
  nbeta = numel(SPM.xX.name);
  %-------%
  
  convec = zeros(ncon, nbases, nbeta);
  
  for c = 1:ncon
    nimg = numel(cfg.evtA(c).img);
    
    for i = 1:nimg
      
      %-get correct weight
      if strfind(cfg.evtA(c).img{i}, 'SlowWave')
        wbeta = wcon{1} / nimg;
      elseif strfind(cfg.evtA(c).img{i}, 'Spindle')
        wbeta = wcon{2} / nimg;
      else
        wbeta = wcon{3} / nimg;
      end
      
      for b = 1:nbases
        %-find beta columns
        if strfind(cfg.evtA(c).img{i}, 'param')
          imgn = [cfg.evtA(c).img{i} '^1']; % add ^1 for parametric
        else
          imgn = [cfg.evtA(c).img{i}];
        end
        
        betaname = sprintf('%s*bf(%1.f)', imgn, b);
        idxCell = strfind(SPM.xX.name, betaname);
        ibeta = ~cellfun(@isempty, idxCell);
        
        %-make vector
        convec(c,b,ibeta) = cfg.evtA(c).con(i) * wbeta;
      end
    end
  end
  %-----------------%
  
  %-----------------%
  %-compute contrasts
  matlabbatch = [];
  
  matlabbatch{1}.spm.stats.con.spmmat = {[gdir 'SPM.mat']};
  matlabbatch{1}.spm.stats.con.delete = 1;
  
  cnt = 0;
  for c = 1:ncon
    for b = 1:nbases
      cnt = cnt + 1;
      cname = sprintf('%s_%1.f', cfg.evtA(c).name, b);
      cvec = squeeze(convec(c, b, :))';

      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = cname;
      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = cvec;
      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.sessrep = 'none';
      
      svec = sprintf(' % 5.3g', cvec);
      outtmp = sprintf('%s: %s\n', cname, svec);
      output = [output outtmp];
    end
  end
  
  spm_jobman('run', matlabbatch)
  %-----------------%
  %---------------------------%
  
  %---------------------------%
  %-copy contrasts
  cnt = 0;
  for c = 1:ncon
    for b = 1:nbases
      cnt = cnt + 1;
      
      for e = 1:2 % img and hdr
        con_from = sprintf('%scon_%04.f.%s', gdir, cnt, ext{e});
        con_to   = sprintf('%s%s_b%02.f_p%02.f.%s', conA, cfg.evtA(c).name, b, subj, ext{e});
        copyfile(con_from, con_to);
      end
    end
  end
  %---------------------------%  
  rmdir(gdir, 's')
  
end
%-------------------------------------%

%-------------------------------------%
%-Specify first-level model: F2B-B2F SW CONTRAST
if strfind(cfg.AorB, 'B')

  %---------------------------%
  %-directories
  %-------%
  %-glm folder
  gdir = [sdir 'glm' filesep];
  
  if exist(gdir, 'dir'); rmdir(gdir, 's'); end
  mkdir(gdir)
  %-------%
  
  %-------%
  conB = [cfg.dirB 'contrasts' filesep];
  
  if ~exist(conB, 'dir'); mkdir(conB); end
  ext = {'hdr', 'img'};
  %-------%
  %---------------------------%
  
  %---------------------------%
  %-specify model
  matlabbatch = [];
  
  %-------%
  %-first level parameters
  matlabbatch{1}.spm.stats.fmri_spec.dir = {gdir};
  matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
  matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.46;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
  %-------%
  
  %-----------------%
  %-loop through recordings (sessions)
  rcnt = 0;
  rall = [];
  
  for r = 1:numel(SW_onset)
    
    if size(SW_onset{r},1) >= cfg.minW && ...
        size(bSW_onset{r},1) >= cfg.minW && ...
        size(sSW_onset{r},1) >= cfg.minW
      rcnt = rcnt + 1; % <- rcnt is the sessions actually written on file (enough SW), r is the sessions from sleep scoring
      rall = [rall r];
      
      %-------%
      %-F2B SW
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).scans = swfIMG{rcnt};
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).name = 'f2bSWmain';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).onset = bSW_onset{r}' + cfg.dlay;
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).duration = bSW_dur{r}';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).tmod = 0;
      if strcmp(cfg.pmod, '')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
      elseif strcmp(cfg.pmod, '_pmod-dur')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).name  = 'f2bSWparam';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).param = bSW_dur{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).poly  = 1;
      elseif strcmp(cfg.pmod, '_pmod-par')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).name  = 'f2bSWparam';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).param = bSW_param{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(1).pmod(1).poly  = 1;
      end
      %-------%
      
      %-------%
      %-B2F SW
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).name = 'b2fSWmain';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).onset = sSW_onset{r}' + cfg.dlay;
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).duration = sSW_dur{r}';
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).tmod = 0;
      if strcmp(cfg.pmod, '')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
      elseif strcmp(cfg.pmod, '_pmod-dur')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).name  = 'b2fSWparam';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).param = sSW_dur{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).poly  = 1;
      elseif strcmp(cfg.pmod, '_pmod-par')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).name  = 'b2fSWparam';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).param = sSW_param{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(2).pmod(1).poly  = 1;
      end
      %-------%
      
      %-------%
      %-Spindles
      if strcmp(cfg.bsSP, '_bsSP')
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(3).name = 'Spindles';
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(3).onset = SP_onset{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(3).duration = SP_dur{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(3).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
      end
      %-------%
      
      if strcmp(cfg.heart, 'yes')
         matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress.name = 'heart';
         matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress.val  = RR{r};
      else
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).regress = struct('name', {}, 'val', {});
      end
      
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).multi = {''};
      if cfg.cvrp
        matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).multi_reg = rp(rcnt);
      end
      matlabbatch{1}.spm.stats.fmri_spec.sess(rcnt).hpf = 128;
      %-------%
    end
  end
  %-----------------%
  
  %-------%
  %-bases
  if strcmp(cfg.bases, 'hrf')
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = cfg.basopt;
  else
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = cfg.basopt(1);
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order  = cfg.basopt(2);
  end
  %-------%
  
  %-------%
  %-first level parameters
  matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
  matlabbatch{1}.spm.stats.fmri_spec.volt = cfg.volt;
  matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
  matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  %-------%
  
  spm_jobman('run', matlabbatch)
  %---------------------------%
  
  %---------------------------%
  %-Estimate first-level model
  matlabbatch = [];
  
  matlabbatch{1}.spm.stats.fmri_est.spmmat = {[gdir 'SPM.mat']};
  matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
  
  spm_jobman('run', matlabbatch)
  %---------------------------%
  
  %---------------------------%
  %-Create contrasts
  load([gdir 'SPM.mat'])
  
  %-----------------%
  %-prepare contrasts
  %-------%
  %-contrast weight (for each session)
  if strcmp(cfg.wcon, 'yes')
    wcon = []; % weight for the contrast
    wcon{1} = cellfun(@numel, bSW_dur);
    wcon{1} = wcon{1}(rall); % it's important to only keep the sessions which are used (have enough events)
    wcon{1} = wcon{1} / sum(wcon{1}); % rescale to one
    
    wcon{2} = cellfun(@numel, sSW_dur);
    wcon{2} = wcon{2}(rall); % it's important to only keep the sessions which are used (have enough events)
    wcon{2} = wcon{2} / sum(wcon{2}); % rescale to one
    
  else
    wcon{1} = ones(1, numel(rall)) / numel(rall);
    wcon{2} = ones(1, numel(rall)) / numel(rall);
  end
  wcon{3} = ones(1, numel(rall)) / numel(rall); % for all the other cases
  %-------%
  
  %-------%
  %-number of bases, contrasts, betas
  if strcmp(cfg.bases, 'hrf');
    nbases = sum(cfg.basopt) + 1;
  else
    nbases =  cfg.basopt(2);
  end
  ncon = numel(cfg.evtB);
  nbeta = numel(SPM.xX.name);
  %-------%
  
  convec = zeros(ncon, nbases, nbeta);
  
  for c = 1:ncon
    nimg = numel(cfg.evtB(c).img);
    
    for i = 1:nimg
      
      %-get correct weight
      if strfind(cfg.evtB(c).img{i}, 'f2b')
        wbeta = wcon{1} / nimg;
      elseif strfind(cfg.evtB(c).img{i}, 'b2f')
        wbeta = wcon{2} / nimg;
      else
        wbeta = wcon{3} / nimg;
      end
      
      for b = 1:nbases
        %-find beta columns
        if strfind(cfg.evtB(c).img{i}, 'param')
          imgn = [cfg.evtB(c).img{i} '^1']; % add ^1 for parametric
        else
          imgn = [cfg.evtB(c).img{i}];
        end
        
        betaname = sprintf('%s*bf(%1.f)', imgn, b);
        idxCell = strfind(SPM.xX.name, betaname);
        ibeta = ~cellfun(@isempty, idxCell);
        
        %-make vector
        convec(c,b,ibeta) = cfg.evtB(c).con(i) * wbeta;
      end
    end
  end
  %-----------------%
  
  %-----------------%
  %-compute contrasts
  matlabbatch = [];
  
  matlabbatch{1}.spm.stats.con.spmmat = {[gdir 'SPM.mat']};
  matlabbatch{1}.spm.stats.con.delete = 1;
  
  cnt = 0;
  for c = 1:ncon
    for b = 1:nbases
      cnt = cnt + 1;
      cname = sprintf('%s_%1.f', cfg.evtB(c).name, b);
      cvec = squeeze(convec(c, b, :))';

      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = cname;
      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = cvec;
      matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.sessrep = 'none';
      
      svec = sprintf(' % 5.3g', cvec);
      outtmp = sprintf('%s: %s\n', cname, svec);
      % output = [output outtmp];
    end
  end
  
  spm_jobman('run', matlabbatch)
  %-----------------%
  %---------------------------%
  
  %---------------------------%
  %-copy contrasts
  cnt = 0;
  for c = 1:ncon
    for b = 1:nbases
      cnt = cnt + 1;
      
      for e = 1:2 % img and hdr
        con_from = sprintf('%scon_%04.f.%s', gdir, cnt, ext{e});
        con_to   = sprintf('%s%s_b%02.f_p%02.f.%s', conB, cfg.evtB(c).name, b, subj, ext{e});
        copyfile(con_from, con_to);
      end
    end
  end
  %---------------------------%
   
  %  rmdir(gdir, 's')
  
end
%-------------------------------------%

%---------------------------%
%-remerge files into 4D structure
for r = 1:numel(swfIMG)
  merge_into_4D(swfIMG{r});
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

function n_vol = count_volumes_from_name(img_name)

firstvol = str2double(img_name(end-12:end-9));
lastvol = str2double(img_name(end-7:end-4));
n_vol = lastvol - firstvol + 1;

function n_vol = count_volumes_in_nii(img_name)

n_vol = bash(['fslinfo ' img_name ' | awk ''NR==5'' | awk ''{print $2}''']);
n_vol = str2double(n_vol);

function output = bash(command, cwd)

if nargin == 1
  [~, output] = system(['. ~/.bashrc; ' command]);
else
  [~, output] = system(['. ~/.bashrc; cd ' cwd ' ; ' command]);
end
