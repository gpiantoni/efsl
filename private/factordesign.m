function factordesign(cfg, conimg, type)
%FACTORDESIGN make quick contrast
% you should use fd (factorial design) instead of anovaw (within-subject
% anova), if you're comparing against baseline. If you want to compare two
% conditions, you compare them at the first-level and then you check if the
% difference is different from zero at the second-level. In theory (with
% linear mixed models, doing the difference at the first or second level
% should not change, as in the case of paired t-test). However, the degrees
% of freedrom are different and the error term combines the intrasubject
% error with the between-subject error (it's still valid, but I think that
% comparing at the first level is better).

mversion = 8;
%08 11/12/13 use the same function for SW and b2fSW
%07 11/12/13 don't use anovaw, bc we're testing each bin against baseline, not against each other, so we don't need to model the subject-specific effect
%06 11/12/12 use simpler full-factorial (like in 
%05 11/12/12 change log (don't use cfg.fid)
%04 11/09/13 catch error for no significant voxels in model estimation
%03 11/09/01 more flexible format
%02 11/04/26 fixed output
%01 11/04/21 created

%---------------------------%
%-directory
%-----------------%
%-SW or F2B-B2F
if type == 1
  dirAB = cfg.dirA;
  spmAB = cfg.spmA;
else
  dirAB = cfg.dirB;
  spmAB = cfg.spmB;
end
%-----------------%

%-------%
%-name is composed of images to contrast
condir = [dirAB conimg.name filesep];
%-------%

if isdir(condir); rmdir(condir, 's'); end
mkdir(condir)

%-----------------%
%-start log
output = sprintf('\n%s (v%02.f) on %s in %s\n', ...
  mfilename,  mversion, conimg.name, spmAB);
%-----------------%
%---------------------------%

%---------------------------%
%-define the contrasts
%-------%
%-find bases
bases = fieldnames(conimg.b);
%-------%

%-------%
%-output
outtmp = sprintf('%s, using %1.f bases\n', ...
  conimg.name, numel(bases));
output = [output outtmp];
%-------%
%---------------------------%

%-------------------------------------%
%-Second-level factorial design
%---------------------------%
matlabbatch = [];
matlabbatch{1}.spm.stats.factorial_design.dir = {condir};

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.name = 'bases';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.levels = numel(bases);
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact.ancova = 0;

for i = 1:numel(bases)
  matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = i;
  matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = conimg.b.(bases{i});
end

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;

%-----------------%
%-define neuroanat mask
switch cfg.msk2
  case {'no' 'con1'}
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
  otherwise
    maskfile = [cfg.mdir cfg.msk2 '.img'];
    if exist(maskfile, 'file')
      matlabbatch{1}.spm.stats.factorial_design.masking.em = {[cfg.mdir cfg.msk2 '.img,1']};
    else
      error(['could not find mask: ' maskfile])
    end
    
    %-------%
    %-output
    ROInames = sprintf(' %s,', cfg.ROIs{:});
    outtmp = sprintf('contrast is masked by %s (%s: %s)\n', cfg.msk2, cfg.atls, ROInames);
    output = [output outtmp];
    %-------%
end
%-----------------%

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch)
%-------------------------------------%

%---------------------------%
%-Second-level factorial estimation
matlabbatch = [];

matlabbatch{1}.spm.stats.fmri_est.spmmat = {[condir 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
try
  spm_jobman('run', matlabbatch)
catch % probably non-significant voxels
  %-------%
  %-output
  output = sprintf(['error in model estimation: probably no significant voxels in ' conimg.name '\n']);
  fprintf(output)
  fid = fopen(cfg.log, 'a');
  fwrite(fid, output);
  fclose(fid);
  %-------%
  return
end
%---------------------------%

%---------------------------%
%-Create contrasts
matlabbatch = [];

matlabbatch{1}.spm.stats.con.spmmat = {[condir 'SPM.mat']};

%-----------------%
%-contrast
conmat = diag(ones(numel(bases),1));
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = conimg.name;
matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {conmat}';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
%-----------------%

matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch)
%---------------------------%

%---------------------------%
%-print results
%-----------------%
%-all the contrasts
matlabbatch = [];

matlabbatch{1}.spm.stats.results.spmmat = {[condir 'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = [spmAB conimg.name];
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none'; % 'FWE' or 'none'
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = true;

spm_jobman('run', matlabbatch)
%-----------------%

%---------------------------%
%-create pdf
pdfname = [datestr(now, 'yymmdd') spmAB '_' conimg.name '_' cfg.msk2 '.pdf'];

system(['ps2pdf spm_' datestr(now, 'yyyymmmdd') '.ps ' pdfname]);

delete(['spm_' datestr(now, 'yyyymmmdd') '.ps'])
copyfile([condir pdfname], [cfg.outp pdfname])
copyfile([condir filesep pdfname], [dirAB pdfname])
%---------------------------%

%--------------------------------------%

%---------------------------%
%-end log
%-----------------%
fprintf(output)
fid = fopen(cfg.log, 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%