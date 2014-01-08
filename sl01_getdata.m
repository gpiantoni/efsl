function sl01_getdata(cfg, subj)
%SLEEPLIEGE 01: getdata

mversion = 11;
%11 11/09/09 2nd argument for subj (and cfg.subj -> subj)
%10 11/08/26 complete data structure, no fmri at all
%09 11/07/28 copy the mat, ln the dat
%08 11/07/28 subjcode is not necessary anymore (convert_data in rec already fixed that)
%07 110705 use dir struct with raw and preproc dir within each subj dir
%06 110705 file name should not change yet
%05 11/05/25 use ln instead of cp
%04 11/03/07 ask user before overwritting folders
%03 11/02/28 uniform cfg
%02 11/02/07 no for subj
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
%-Copy data from raw
disp(subj)
recd = sprintf('%s%04.f/', cfg.recs, subj);
subd = sprintf('%s%04.f/', cfg.data, subj);
if ~isdir(subd); mkdir(subd); end

%-------%
%-rawd and prep for each subj
% subjcode = [809 845 862 944 990 1153 1156 1304 1309 1337 1351 1380 1382 1402];
%-------%

if strcmpi(cfg.type, 'fmri') || strcmpi(cfg.type, 'both')
  disp('fmri')
  
% This probably does not work anymore. Use: sl05_divide_rec to get hdr/img files  
%   %-------%
%   %-fmri
%   ana = [subd 'fmri' filesep cfg.rawd]; % from rawdata
%   fmri = [cfg.data ndir filesep 'fmri' filesep];
%   if exist(fmri, 'dir')
%     ButtonName = questdlg([fmri ' already exists, do you want to overwrite it?'], ...
%                          'Overwrite fmri folder', ...
%                          'Yes', 'No', 'No');
%     if ~strcmp(ButtonName, 'Yes')
%       warning(['using old folder: ' fmri])
%       return
%     end
%     rmdir(fmri, 's')
%   end
% 
%   mkdir(fmri)
%   
%   % copy only images and hdr (not the mat)
%   % use bash bc copyfile copies one file at the time
%   system(['ln ' ana 'f*' fmri]);
%   %-------%
% 
%   %-------%
%   %-mri
%   ana = [subd 'mri' filesep cfg.rawd]; % from rawdata
%   mri = [cfg.data ndir filesep 'mri' filesep];
% 
%   rmdir(mri, 's')
%   mkdir(mri)
%   
%   % copy structural (3 localizer, 1 structural)
%   % but not the swf
%   system(['ln ' ana 's*.img ' mri]);
%   %-------%
  
end

if strcmpi(cfg.type, 'eeg') || strcmpi(cfg.type, 'both')
  disp('eeg')
  inp = [recd cfg.eegd filesep cfg.rawd]; % from rawdata
  eeg = [subd cfg.eegd filesep];
  
  if exist(eeg, 'dir')
    ButtonName = questdlg([eeg ' already exists, do you want to overwrite it?'], ...
                         'Overwrite EEG folder', ...
                         'Yes', 'No', 'No');
    if ~strcmp(ButtonName, 'Yes')
      warning(['using old folder: ' eeg])
      return
    end
    rmdir(eeg, 's')
  end
  mkdir(eeg)
  
  %-------%
  system(['ln -s ' inp '*.dat ' eeg]);
  system(['cp ' inp '*.mat ' eeg]);
  system(['chmod u+w ' eeg '*.mat']);
  %-------%
  
  %-------%
  %-fix name in FASST
  allmat = dir([eeg '*.mat']);
  for i = 1:numel(allmat)
    load([eeg allmat(i).name])
    D.fname = allmat(i).name;
    D.path = eeg;
    D.data.fnamedat = [allmat(i).name(1:end-3) 'dat'];
    D.data.y.fname = [eeg allmat(i).name(1:end-3) 'dat'];
    D = meeg(D);
    save(D);
    clear D
  end
  %-------%

end

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
