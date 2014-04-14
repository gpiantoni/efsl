function sl13_report(cfg)
%SLEEPLIEGE 13: REPORT INTERESTING INFO
% -activation in SW: map, p-value
% -activation in f2b SW: map, p-value, estimated effect
% -activation in b2f SW: map, p-value, estimated effect
% -estimated effect in rLC and lLC
% -activation in PPI: map, p-value

mversion = 7;
%07 12/01/06 write to csv even if not significant
%06 11/12/23 write to csv
%05 11/12/21 report for sw with the "everything-in-contrast" format
%04 11/09/22 include minimum amount of traveling for slow waves
%03 11/09/13 use 'ResMS.hdr' instead of SPM.mat (it checks if the model was estimated)
%02 11/09/02 report results as in spm_list, with ROI and BA
%01 11/09/01 created

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

output = sprintf('%s\n\n=======================================\nREPORT\n', output);

%---------------------------%
%recap of slow waves
output = sprintf('%s===================\n', output);
output = sprintf('%sSlow wave characteristics\n', output);

%-----------------%
%-loop over subjects
for s = 1:14
  trdir  = sprintf('%s%04.f%s', cfg.data, s, '/spm/triggers/');
  output = sprintf('%sp%02.f:\n', output, s);
  
  %--------%
  %-load triggers
  load([trdir cfg.trigA], 'SW_onset')
  load([trdir cfg.trigB], 'bSW_param', 'sSW_param')
  %--------%
  
  %--------%
  %-loop over sessions
  for ss = 1:numel(SW_onset)
    output = sprintf('%s   in s%1.f: nSW% 4.f', output, ss, numel(SW_onset{ss}));
    output = sprintf('%s (f2b:% 4.f, % 4.f;', output, numel(bSW_param{ss}), mean(bSW_param{ss}));
    output = sprintf('%s b2f:% 4.f, % 4.f)\n', output, numel(sSW_param{ss}), mean(sSW_param{ss}));
  end
  %--------%
  
end
%-----------------%

output = sprintf('%s===================\n', output);
%---------------------------%


%---------------------------%
%recap of slow waves
output = sprintf('%s===================\n', output);
output = sprintf('%sSlow wave summary\n', output);

output = [output sw_summary(cfg)];

output = sprintf('%s===================\n', output);
%---------------------------%

%---------------------------%
%activation in main SW: map, p-value, ROI
for i = 1:numel(cfg.evtA)
  dirAcon = [cfg.dirA cfg.evtA(i).name filesep];
  
  if exist([dirAcon 'ResMS.hdr'], 'file')
    output = sprintf('%s===================\n', output);
    output = sprintf('%s%s\n', output, cfg.evtA(i).name);
    output = [output peak_n_ROI(dirAcon, 'FWE',  1, cfg.gdLC)];
    output = [output peak_n_ROI(dirAcon, 'none',  1)];
    output = sprintf('%s===================\n\n', output);
  end
end
%---------------------------%

%---------------------------%
%-write to csv
fid = fopen(cfg.csvf, 'a+');
fprintf(fid, struct2log(cfg, 'csv'));
%---------------------------%

%---------------------------%
%activation in b2fSW and f2bSW: map, p-value, ROI
for i = 1:numel(cfg.evtB)
  dirBcon = [cfg.dirB cfg.evtB(i).name filesep];
  
  if exist([dirBcon 'ResMS.hdr'], 'file')
    output = sprintf('%s===================\n', output);
    output = sprintf('%s%s\n', output, cfg.evtB(i).name);
    [outtmp outcsv] = peak_n_ROI(dirBcon, 'FWE',  1, cfg.gdLC);
    output = [output outtmp];
    output = [output peak_n_ROI(dirBcon, 'none',  1)];
    output = sprintf('%s===================\n\n', output);
  else
    outcsv = ['no,no,'];
  end
  fprintf(fid, outcsv);
end
%---------------------------%

%---------------------------%
%activation in PPI: map, p-value
dirPPI = [cfg.dirB cfg.evtB(cfg.LCic).name filesep 'PPI_' cfg.LCdf filesep];
if isdir(dirPPI)
  output = sprintf('%s===================\n', output);
  output = sprintf('%sPPI SW\n', output);
  output = [output peak_n_ROI(dirPPI, 'none', 5)];
  output = sprintf('%s===================\n\n', output);
end
%---------------------------%

%-------%
output = sprintf('%s=======================================\n\n', output);
fprintf(output)

fprintf(fid, sprintf('\n'));
fclose(fid);
%-------%

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