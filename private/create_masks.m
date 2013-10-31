function create_masks(cfg)
%CREATE_MASKS
% create_masks(cfg, ROI, cfg.msk2)
%
% cfg
%  .Dwfu = dir of WFU templates
%  .dMsk = dir of masks
%
%  .atls = 'aal_MNI_V4'
%  .ROIs = ROIs (has to be a cell)
%
%  .msk2 = file name it saves the mask in.

% 11/09/02 added masksize
% 11/04/27 allows for hand-drawn ROI
% 11/02/10 created

%-----------------%
%-does it want a mask at all?
output = '';
if any(strcmp(cfg.msk2, {'no', 'con1'}))
  return
  
elseif ~isempty(strfind(cfg.msk2, 'drawn'))
  
  %-------%
  %-output
  outtmp = sprintf('drawn mask %s has size: %5.2fcm3\n', cfg.msk2, masksize([cfg.dMsk cfg.msk2 '.img']));
  output = [output outtmp];
  %-------%
  
  %-----------------%
  fprintf(output)
  fid = fopen(cfg.log, 'a');
  fwrite(fid, output);
  fclose(fid);
  %-----------------%
  
  return
  
end
%-----------------%

%---------------------------%
%-input check
addpath('/data/toolbox/WFU_PickAtlas_3.0.1/wfu_pickatlas/')

descrip = ['mask ' cfg.msk2 ': ' sprintf('%s,', cfg.ROIs{:}) ' in ' cfg.atls];

%-----------------%
%- if it exists, check that the name matches the content
if exist([cfg.dMsk cfg.msk2 '.img'], 'file')
  hMsk = spm_vol([cfg.dMsk cfg.msk2 '.img']); % header
  
  if strcmp(hMsk.descrip, descrip)
    %-------%
    %-output
    outtmp = sprintf('mask: %s already exists and matches descrip (size: %5.2fcm3)\n\n', cfg.msk2, masksize([cfg.dMsk cfg.msk2 '.img']));
    output = [output outtmp];
    %-------%
    return
  else
    %-------%
    %-output
    outtmp = sprintf('mask: %s already exists but doesn''t match descrip.\nOld descrip: %s\nNew descrip: %s\nIt will be over-written\n', ...
      cfg.msk2, hMsk.descrip, descrip);
    output = [output outtmp];
    %-------%
  end
  
end
%-----------------%
%---------------------------%

%---------------------------%
%-read vol and txt
hAtl = spm_vol([cfg.Dwfu cfg.atls '.nii']); % header
vAtl = spm_read_vols(hAtl); % volume

allROI = wfu_txt2roi( [cfg.Dwfu cfg.atls '.txt']); % <- part of WFU toolbox

iAtl = [allROI.ID]; % index of Atlas
lAtl = {allROI.Nom_C}; % labels of Atlas
%---------------------------%

%---------------------------%
%-prepare mask image

hMsk = hAtl;
vMsk = zeros(hAtl.dim);

hMsk.fname = [cfg.dMsk cfg.msk2 '.img'];
hMsk.descrip = descrip;

%-----------------%
%-add ROI
for r = 1:numel(cfg.ROIs)
  iLbl = strcmpi(lAtl, cfg.ROIs{r});  % index of labels
  
  if isempty(iLbl) || numel(find(iLbl)) > 1
    error(['ROI name ' cfg.ROIs{r} ' doesn''t match labels'])
  end
  
  iROI = iAtl(iLbl); % they don't always correspond
  
  vMsk = vMsk + (vAtl == iROI);
  
end

spm_write_vol(hMsk, vMsk);
%-----------------%

%-------%
%-output
outtmp = sprintf('created %s (size: %5.2fcm3)\n\n', descrip, masksize([cfg.dMsk cfg.msk2 '.img']));
output = [output outtmp];
%-------%

%-----------------%
fprintf(output)
fid = fopen(cfg.log, 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%

function cubcm = masksize(img)

hmask = spm_vol(img); % header
vmask = spm_read_vols(hmask);

cubcm = numel(find(vmask(:))) * abs(det(hmask.mat)) / 1e3; % from mm3 -> cm3

