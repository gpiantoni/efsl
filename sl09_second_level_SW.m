function sl09_second_level_SW(cfg)
%SLEEPLIEGE 09: SECOND LEVEL SW type

mversion = 12;
%12 11/12/12 collect the contrasts in different ways
%11 11/09/01 cfg.evtB is a struct, 'collect the con images' only collects images, the contrasts are run within factordesign (more flexible)
%10 11/07/15 cfg.evtB is in final_efsl.m (it can be reused)
%09 11/04/26 simplified collection of contrasts
%08 11/04/26 uses peak_in_LC to check whether the peak is in LC
%07 11/04/26 uses factor design (no useless columns in design matrix)
%06 11/04/13 convert -> ps2pdf
%05 11/04/01 depreceated cfg.mask2 (never use contrast 1 for the other contrasts)
%04 11/02/10 only for contrast big/small slow waves
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
if ~exist(cfg.dirB, 'dir')
  error(['directory ' cfg.dirB ' does not exist'])
end

conB = [cfg.dirB 'contrasts' filesep];
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
for k = 1:numel(cfg.evtB) % contrast groups
  for b = 1:nbases % bases of each image
    bname = sprintf('b%02.f', b);
    
    conpat = sprintf('%s%s_%s_p*.img', ...
      conB, cfg.evtB(k).name, bname); % contrast pattern
    
    allcon = dir(conpat);
    
    for s = 1:numel(allcon)
      conIMG(k).name = cfg.evtB(k).name;
      conIMG(k).b.(bname){s,1}= [conB allcon(s).name ',1'];
    end
  end
end
%-----------------%
%---------------------------%

%---------------------------%
%-FACTOR DESIGN
for i = 1:numel(conIMG)
   factordesign(cfg, conIMG(i), 2)
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
