function sl08_second_level_main(cfg)
%SLEEPLIEGE 08: SECOND LEVEL MAIN

mversion = 9;
%10 11/12/13 use factoranalysis
%09 11/08/31 contrasts are specified in main file, run only with one level
%08 11/08/31 use mask2 for the main contrast as well
%07 11/04/13 convert -> ps2pdf
%06 11/03/29 fix contrasts, allows for parametric (specified) and spindles (not specified) contrasts
%05 11/02/16 allows for FIR, flexible contrasts. But it's more tricky to collect the contrasts
%04 11/02/10 only for main analysis
%03 11/01/27 print output, into pdf and send it to general pdf
%02 11/01/27 removed try, add log
%01 11/01/18 created

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-Collect contrasts
%-----------------%
%-prepare directories
if ~exist(cfg.dirA, 'dir')
  error(['directory ' cfg.dirA ' does not exist'])
end

conA = [cfg.dirA 'contrasts' filesep];
%-----------------%

%-----------------%
%-names
if strcmp(cfg.bases, 'hrf');
  nbases = sum(cfg.basopt) + 1;
else
  nbases =  cfg.basopt(2);
end
%-----------------%

%-----------------%
%-collect the con images
conIMG = [];
for k = 1:numel(cfg.evtA) % contrast groups
  for b = 1:nbases % bases of each image
    bname = sprintf('b%02.f', b);
    
    conpat = sprintf('%s%s_%s_p*.img', ...
      conA, cfg.evtA(k).name, bname); % contrast pattern
    
    allcon = dir(conpat);
    
    for s = 1:numel(allcon)
      conIMG(k).name = cfg.evtA(k).name;
      conIMG(k).b.(bname){s,1}= [conA allcon(s).name ',1'];
    end
  end
end
%-----------------%
%---------------------------%

%---------------------------%
%-FACTOR DESIGN
for i = 1:numel(conIMG)
   factordesign(cfg, conIMG(i), 1)
end
%---------------------------%

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
