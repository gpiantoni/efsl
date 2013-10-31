function sl11_ppisubj(cfg, subj)
%SLEEPLIEGE 11: PPISUBJ
% Do GLM on PPI design matrix for each subject
% My version (in /scripts/progress/) was more flexible bc I tried to take
% care of: 
%   * neural activity in LC should be estimated not using HRF (requires
%     large changes in the spm8 code)
%   * It should be possible to use hrf + derivative for the PPI, instead of
%     only one column (partially implemented in my code)

mversion = 4;
%04 11/09/23 2nd argument for subj (and cfg.subj -> subj)
%03 11/09/01 using cfg.evtB(cfg.LCic).name instead of cfg.cdir
%02 11/07/18 added filesep to cfg.cdir
%01 11/07/15 using PPPI

spm_jobman('initcfg')

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

%-----------------%
%-directories
dirP = [cfg.dirB cfg.evtB(cfg.LCic).name filesep 'PPI_' cfg.LCdf filesep];
conP = [dirP 'contrasts' filesep];
if ~exist(conP, 'dir'); mkdir(conP); end
  
%-------%
%-subject
ndir = sprintf('%04.f', subj);
sdir = sprintf('%s%s/spm/', cfg.data, ndir);
gdir = [sdir 'glm' filesep];
pdir = [gdir 'PPI_' cfg.LCdf filesep];

if isdir(pdir); rmdir(pdir, 's'); end
%-------%
%-----------------%

%-------------------------------------%
%---------------------------%
%-PPPI
P = [];

%-----------------%
%-parameters
P.Region = cfg.LCdf;

%-------%
%-ROI
if strcmp(cfg.LCdf(1:3), 'ROI') 
  P.VOI = [cfg.mdir cfg.LCdf];
else
  P.VOI = [cfg.dirB cfg.evtB(cfg.LCic).name filesep 'ROI_' cfg.LCdf '.img'];
end
%-------%

P.subject = ndir;
P.directory = gdir;

P.extract = cfg.pext; 
P.method = cfg.pmet;
P.Tasks = cfg.evtB(cfg.ptsk).img;
%-----------------%

%-----------------%
%-default
P.Estimate = 1;
P.contrast = {'Omnibus F-test for PPI Analyses'};
P.Weights = ones(1,numel(P.Tasks));
% if strcmp(cfg.pmet, 'cond'); P.Weights = ones(numel(P.Tasks)/2, 1); end
P.maskdir = [];
P.equalroi = 1;
P.FLmask = 0;
P.VOI2 = {};
P.analysis = 'psy';
P.SPMver = 8;
P.CompContrasts = 0; % I'll compute the contrasts myself
P.Weighted = [];

P.Contrasts = [];
%-----------------%

PPPI(P);
%---------------------------%

%---------------------------%
%-compute contrasts and move them
%-----------------%
%-------%
SPMfile = [pdir 'SPM.mat'];
load(SPMfile, 'SPM')

convec = ~cellfun(@isempty, strfind({SPM.Vbeta.descrip}, 'PPI'));
%-------%

matlabbatch = [];

matlabbatch{1}.spm.stats.con.spmmat = {SPMfile};
matlabbatch{1}.spm.stats.con.delete = 1;

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'PPI';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = double(convec);
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

spm_jobman('run', matlabbatch)
%-----------------%

%-----------------%
%-Copy contrasts
ext = {'hdr', 'img'};

cname = ['con_ppi' cfg.evtB(cfg.ptsk).name];

for e = 1:2 % img and hdr
  con_from = sprintf('%scon_%04.f.%s', pdir, 1, ext{e});
  con_to   = sprintf('%s%s_p%02.f.%s', conP, cname, subj, ext{e});
  copyfile(con_from, con_to);
end
%-----------------%
%---------------------------%
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