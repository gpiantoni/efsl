function sl10_roi4ppi(cfg)
%SLEEPLIEGE 10: ROI4PPI
% define extra column for PPI (LC and SW)

mversion = 6;
%06 11/12/13 onlyLC: only select clusters that include LC
%05 11/09/01 using cfg.evtB(cfg.LCic).name instead of cfg.cdir
%04 11/07/20 allows for biact
%03 11/07/18 added filesep to cfg.cdir
%02 11/07/15 rewrite using donald mcdonald peak_nii
%01 11/03/29 new implementation, based on spm_regions

% cfg.LCdf 'act' (activation, FWE, .05) or 'biact' (make act bilateral)
% maxclu (only largest cluster), bimaxclu (make maxclu bilateral)

%---------------------------%
%-start log
output = sprintf('%s (v%02.f) started at %s on %s\n', ...
  mfilename,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%-------------------------------------%
%-definition of ROI (LC)
if strcmp(cfg.LCdf, 'act') || strcmp(cfg.LCdf, 'biact')
  
  %---------------------------%
  %-get peak activation
  xSPM = [];
  
  xSPM.swd       = [cfg.dirB cfg.evtB(cfg.LCic).name filesep];
  xSPM.Ic        = 1; % index of the contrasts (SW type or bSW)
  xSPM.u         = cfg.LCpv;
  xSPM.Im        = [];
  xSPM.thresDesc = cfg.LCcr; 
  xSPM.title     = 'SW_type';
  xSPM.k         = 0;
  xSPM.units = {'mm' 'mm' 'mm'};
  
  [SPM, xSPM] = spm_getSPM(xSPM);
  
  %-----------------%
  %-only use cluster within the LC
  if strcmp(cfg.onlyLC, 'yes')
    
    [~, iLC] = intersect(xSPM.XYZmm', cfg.gdLC, 'rows');
    spmcl = spm_clusters(xSPM.XYZ); % define clusters
    
    XYZ = []; XYZmm = [];
    for i = 1:numel(iLC)
      inLC = spmcl == spmcl(iLC(i)); % this will break if two clusters are active
      XYZ = [XYZ xSPM.XYZ(:,inLC)];
      XYZmm = [XYZmm xSPM.XYZmm(:,inLC)];
    end
    
  else
    XYZ = xSPM.XYZ;
    XYZmm = xSPM.XYZmm;
  end
  %-----------------%
  
  %-----------------%
  %-abort if no voxels are significant
  if size(XYZ, 2) == 0
    
    %-------%
    %-output
    outtmp = sprintf('there are no peak voxel, aborting\n');
    output = [output outtmp];
    %-------%
    
    return
  end
  %-----------------%
  %---------------------------%
  
  %---------------------------%
  % other options to specify LC
  if strcmp(cfg.LCdf, 'biact')
    XYZbi = XYZmm .* repmat([-1; 1; 1], 1, size(XYZmm, 2));
    XYZmm = [XYZmm XYZbi];
    
    XYZ = SPM.xVol.iM * [XYZmm; ones(1, size(XYZmm,2))];
    XYZ(4,:) = [];
  end
  
  if strcmp(cfg.LCdf, 'maxclu') % largest cluster
    %-identify cluster with peak voxel
    error('not implemented yet')
    
    % [~, imaxZ] = max(xSPM.Z);
    % spmcl = spm_clusters(xSPM.XYZ); % define clusters
    % icl = spmcl(imaxZ);
    %
    % cl_mm  = xSPM.XYZmm(:, spmcl == icl);
    % cl_i = xSPM.XYZ(:, spmcl == icl);
  end
  
  if strcmp(cfg.LCdf, 'bimaxclu')
    error('not implemented yet')
  end
  %---------------------------%
  
  %---------------------------%
  %-report voxels
  for v = 1:size(XYZ, 2)
    %-------%
    %-output
    outtmp = sprintf('mask includes [%1.f, %1.f, %1.f]\n', ...
      XYZmm(1, v), XYZmm(2, v), XYZmm(3, v));
    output = [output outtmp];
    %-------%
  end
  %---------------------------%
  
end
%-------------------------------------%

%-------------------------------------%
%-save to file
%---------------------------%
%-write ROI to file
ROIn = ['ROI_' cfg.LCdf];

%-------%
%-based on mask
V = spm_vol([cfg.dirB cfg.evtB(cfg.LCic).name filesep 'mask.hdr']);
%-------%

%-------%
Y = zeros(V.dim);
for v = 1:size(XYZ,2)
  Y(XYZ(1,v), XYZ(2,v), XYZ(3,v)) = 1;
end
%-------%

%-------%
%-write to file
V.fname = [ROIn '.img'];
V.private.dat.fname = V.fname;
spm_write_vol(V, Y);
%-------%
%---------------------------%

%---------------------------%
%-save to mat
save([cfg.dirB cfg.evtB(cfg.LCic).name filesep ROIn '.mat'], 'XYZ', 'XYZmm')
%---------------------------%
%-------------------------------------%

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