function sl05_divide_rec(cfg, subj)
%SLEEPLIEGE 05: DIVIDE_REC
%  sl01_getdata might not be necessary for fmri (only eeg)

mversion = 10;
%10 11/09/22 include minimum amount of traveling for slow waves
%09 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%08 11/07/05 get data from recordings and put them in projects
%07 11/02/18 names change depending on smoothing (although smoothing is applied later)
%06 11/02/11 keep all struct together, in cfg.wmsi, not in rdir
%05 11/02/10 include minW (min amount of slow wave to call it a session)
%04 11/02/08 uniform cfg
%03 11/02/07 no for subj
%02 11/01/16 works on SomerenServer
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

%-----------------%
%-files and dir
sIMG = sprintf('%ss-%s_%04.f.img', cfg.wmsd, cfg.code, subj);
mdir = sprintf('%s%04.f/smri/%s', cfg.recs, subj, cfg.rawd);
fdir = sprintf('%s%04.f/fmri/%s', cfg.recs, subj, cfg.rawd);
rdir = sprintf('%s%04.f%s', cfg.data, subj, '/rec/');
if isdir(rdir); rmdir(rdir, 's'); end
mkdir(rdir)
ext = {'img' 'hdr'};
load(cfg.mrkr, 'mkr')

% create sessinfo
for s = [1 3:8 10:12]
  sessinfo(s).sess = 3;
end
for s = [9 13]
  sessinfo(s).sess = 4;
end
%-----------------%

%-----------------%
if ~exist(sIMG, 'file')
  
  %-------%
  %-output
  outtmp = sprintf('%s doesn''t exist, moving structural to %s\n', ...
    sIMG, cfg.wmsd);
  output = [output outtmp];
  %-------%
  
  slist = dir(sprintf('%s%s_%04.f_smri_sleep*img', mdir, cfg.code, subj));
  
  % (p01 has low-res str, while p05 has very high-res struct)
  [dum, smax] = max( [slist.bytes] ); % the largest s is the main structural
  for e = 1 : 2 % img and hdr extensions
    file1 = [mdir slist(smax).name(1:end-3) ext{e}];
    file2 = sprintf('%ss-%s_%04.f.%s', cfg.wmsd, cfg.code, subj, ext{e});
    system(['cp ' file1 ' ' file2]);
    system(['chmod u+w ' file2]);
  end
  
else
  
  %-------%
  %-output
  outtmp = sprintf('%s already exists, not copying structural\n', ...
    sIMG);
  output = [output outtmp];
  %-------%
  
end
%-----------------%

%-----------------%
%-load triggers to see if there are enough slow waves
trdir  = sprintf('%s%04.f%s', cfg.data, subj, '/spm/triggers/');

load([trdir cfg.trigA], 'SW_onset')
load([trdir cfg.trigB], 'bSW_onset', 'sSW_onset')
%-----------------%

%-----------------%
%-functional
rcnt = 0;
for r = 1 : size( mkr(subj).mkr, 1)
  fprintf('p%02.f %02.f/%02.f\n', subj, r, size( mkr(subj).mkr, 1));
  
  %-------%
  %-don't use sessions if too few scans
  if size(SW_onset{r},1) >= cfg.minW && ...
      size(bSW_onset{r},1) >= cfg.minW && ...
      size(sSW_onset{r},1) >= cfg.minW
    rcnt = rcnt + 1; % <- rcnt is the sessions actually written on file (enough SW), r is the sessions from sleep scoring
    
    sdir = sprintf('%sr%02.f/', rdir, rcnt);
    mkdir(sdir)
    
    scan = mkr(subj).mkr(r,1):mkr(subj).mkr(r,2);
    
    for k = 1 : numel(scan)
      oldscan = scan(k);
      
      if subj == 2
        if oldscan < 535 % in mkr, it's always like that, but for consistency
          sesssubj = 3;
        else
          sesssubj = 4;
        end
      elseif subj == 14
        if oldscan < 2629
          sesssubj = 3;
        else
          oldscan = scan(k) - 2629; % to offset the changes in mkr
          sesssubj = 5;
        end
      else
        sesssubj = sessinfo(subj).sess;
      end
      
      for e = 1 : 2 % img and hdr extensions
        file1 = sprintf('%s%s_%04.f_fmri_sleep_s%01.fsc%04.f.%s', ...
          fdir, cfg.code, subj, sesssubj, oldscan, ext{e});
        file2 = sprintf('%sf%02.f-r%02.f-s%1.f-%04.f.%s', ...
          sdir, subj, rcnt, cfg.smoo, scan(k), ext{e});
        
        if strcmp(ext{e}, 'hdr')
          copyfile(file1, file2);
          system(['chmod u+w ' file2]);
        else
          copyfile(file1, file2);
        end
      end
      
    end
    
  else
    %-------%
    %-output
    outtmp = sprintf('r%1.f was skipped (min n SW: %1.f, actual n SW: %1.f)\nThe number of r will not match, no problem\n', ...
      r, cfg.minW, size(SW_onset{r},1));
    output = [output outtmp];
    %-------%
  end
  %-------%
  
end
%-----------------%

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
