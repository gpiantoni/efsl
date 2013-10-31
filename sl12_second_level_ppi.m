function sl12_second_level_ppi(cfg)
%SLEEPLIEGE 12: SECOND LEVEL PPI

mversion = 5;
%05 11/09/01 using cfg.evtB(cfg.LCic).name instead of cfg.cdir
%04 11/07/18 added filesep to cfg.cdir
%03 11/07/15 simpler version based on PPPI (no hrf deriv), old version in scripts/progress
%02 11/04/13 convert -> ps2pdf
%01 11/03/30 created

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-Collect contrasts

%-------%
dirP = [cfg.dirB cfg.evtB(cfg.LCic).name filesep 'PPI_' cfg.LCdf filesep];
%-------%

%-----------------%
%-contrasts name
%-------%
%-prepare directories
conP = [dirP 'contrasts' filesep];
%-------%

cname = ['con_ppi' cfg.evtB(cfg.ptsk).name];

%-----------------%
conIMG = [];

conpat = sprintf('%s%s_p*.img', conP, cname); % contrast pattern;
allcon = dir(conpat);

for s = 1:numel(allcon)
  conIMG.(cname){s,1} = [conP allcon(s).name ',1'];
end
%-----------------%
%---------------------------%

%---------------------------%
%-Second-level factorial design

matlabbatch = [];

matlabbatch{1}.spm.stats.factorial_design.dir = {dirP};

%-----------------%
%- bases factor
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'bases';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = numel(conAll);
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
% %-----------------%
% 
% cnt = 0;
% for k = 1:numel(cond)
%   for c = 1:numel(conAll)
%     cnt = cnt + 1;
%     cname = [cond{k} conAll{c}]; % <- not a cell in this case, but it's very useful to select the right conIMG
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cnt).levels = [c; k];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cnt).scans = conIMG.(cname);
%   end
% end

%-----------------%
% t1 test
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conIMG.(cname);

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run', matlabbatch)
%---------------------------%

%---------------------------%
%-Second-level factorial estimation

matlabbatch = [];

matlabbatch{1}.spm.stats.fmri_est.spmmat = {[dirP 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch)
%---------------------------%

%---------------------------%
%-Create contrasts

matlabbatch = [];

matlabbatch{1}.spm.stats.con.spmmat = {[dirP 'SPM.mat']};

matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'PPI';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {1}';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch)
%---------------------------%

%---------------------------%
%-print results

matlabbatch = [];

matlabbatch{1}.spm.stats.results.spmmat = {[dirP 'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = [cfg.spmB 'ppi'];
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = true;

spm_jobman('run', matlabbatch)
%---------------------------%

%---------------------------%
%-create pdf
pdfname = [datestr(now, 'yymmdd') cfg.spmB '_ppi' '.pdf'];

system(['ps2pdf spm_' datestr(now, 'yyyymmmdd') '.ps ' pdfname]);
delete(['spm_' datestr(now, 'yyyymmmdd') '.ps'])
copyfile([dirP pdfname], [cfg.outp pdfname])
copyfile([dirP pdfname], [cfg.dirB pdfname])
%---------------------------%
%--------------------------------------%

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('%s (v%02.f) ended at %s on %s after %s\n\n', ...
  mfilename, mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen(cfg.log, 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
