function [output outcsv] = peak_n_ROI(condir, corrp, k, ROI)
%PEAK_N_ROI: check whether peak voxel is in ROI and its p-value
%use as output = peak_n_ROI(condir, corrp, k, ROI)
% where
%   condir is a directory with SPM.mat
%   corrp is the type of correction 'none' of 'FWE'
%   k is extent threshold
%   ROI is an MNI point (f.e. LC)

%09 12/01/05 fix outcsv, reports number of clusters and if any of those is in LC
%08 11/12/23 second output to write2csv
%07 11/09/22 compare f2b and b2f contrasts
%06 11/09/02 input needs k instead of Ic (index of the contrast)
%05 11/09/01 renamed peak_n_ROI: tells peaks and distance from ROI (f.e.LC)
%04 11/04/26 added reportpval and report peak voxel of clusters which are more significant than cluster with LC (complicated code)
%03 11/04/26 output checks if LC is significant and all the clusters
%02 11/04/26 requires Ic (contrast index) as input
%01 11/04/26 created

%---------------------------%
%-check input
if ~isdir(condir) || ~exist([condir 'SPM.mat'], 'file')
  error([condir ' is not a directory or SPM.mat does not exist'])
end

if strcmp(corrp, 'none')
  pval = 0.001;
elseif strcmp(corrp, 'FWE')
  pval = 0.05;
end
%---------------------------%

%-------------------------------------%
%-cluster stats
xSPM = [];

xSPM.swd       = condir;
xSPM.Ic        = 1;
xSPM.u         = pval;
xSPM.Im        = [];
xSPM.thresDesc = corrp;
xSPM.title     = condir;
xSPM.k         = k;
xSPM.units = {'mm' 'mm' 'mm'};

[SPM, xSPM] = spm_getSPM(xSPM);
%-------------------------------------%

%-------------------------------------%
%-report results
%-----------------%
%-output
output = sprintf('%s\n%s, extent: %g\n', condir, xSPM.thresDesc, k);
outcsv = '';
%-----------------%

if size(xSPM.Z, 2) == 0
  
  %---------------------------%
  %-no significant voxels
  
  %-----------------%
  %-output
  output = sprintf('%sno significant voxel at %s\n', ...
    output, xSPM.thresDesc);
  outcsv = '0,';
  %-----------------%
  
  %-----------------%
  %-check anyway it is significant in ROI
  if nargin == 4
    outtmp = compare_con(SPM, ROI);
    output = [output outtmp];
    outcsv = sprintf('%s%1.f,', outcsv, numel(strfind(outtmp, sprintf('\n'))));
  end
  %-----------------%
  %---------------------------%
  
else
  
  %---------------------------%
  %-work with clusters
  outcsv = sprintf('%1.f,', size(xSPM.XYZ,2));
  
  %-----------------%
  %-get cluster and sorted Z-values
  cl = spm_clusters(xSPM.XYZ); % define clusters
  [~, Zsort] = sort(xSPM.Z, 'descend');
  clsort = cl(Zsort);
  %-----------------%
  
  %-----------------%
  %-loop over clusters
  %collect info as in the statistics part of SPM (see spm_list)
  %it looks complicated but it should be identical to "results" from SPM
  [~, ncl] = unique(clsort, 'first');
  ncl = clsort(ncl);
  cls = [];
  
  for c = 1:numel(ncl)
    i = ncl(c);
    cls(c).nvox = numel(find(cl==i));
    cls(c).XYZmm = xSPM.XYZmm(:, cl==i)';
    cls(c).XYZ = xSPM.XYZ(:, cl==i)';
    
    peakvox = Zsort(find(clsort == i, 1));
    cls(c).peak = xSPM.XYZmm(:, peakvox)';
    cls(c).Pz = spm_P(1, 0, xSPM.Z(peakvox), xSPM.df, xSPM.STAT,      1, xSPM.n, xSPM.S);
    cls(c).Pu = spm_P(1, 0, xSPM.Z(peakvox), xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.S);
  end
  %-----------------%
  
  %-----------------%
  %-report clusters
  inLC = false;
  for c = 1:numel(cls)
    output = sprintf('%sFWE: %0.3f, uncorr: %0.3f, eff:% 7.2f, nvox:% 3.f, peak [% 4.f % 4.f % 4.f]', ...
      output, cls(c).Pu, cls(c).Pz, 1e3*get_beta(SPM, cls(c).XYZ) ,cls(c).nvox, cls(c).peak(1), cls(c).peak(2), cls(c).peak(3));
    
    if nargin == 4
      %---------%
      %-user asks for ROI
      inROI = intersect(cls(c).XYZmm, ROI, 'rows');
      if ~isempty(inROI) % and it matches!
        output = sprintf('%s  <-- includes ROI [% 4.f % 4.f % 4.f]\n', output, inROI(1,1), inROI(1,2), inROI(1,3)); % only show first one if there are more
        inLC = true;
      else
        output = sprintf('%s\n', output);
      end
      %---------%
      
    else
      %---------%
      %-no ROI: report area of main activation
      atlas = mni2ba(cls(c).XYZmm, 6);
      atlas = count_unique({atlas.aal});
      
      output = sprintf('%s  (%s)\n', output, atlas{1});
      %---------%
    end
  end
  
  %---------%
  if nargin == 4
    if inLC % it at least one cluster is in LC
      outcsv = [outcsv 'inLC,'];
    else
      outcsv = [outcsv 'noLC,'];
    end
  end
  %---------%
  %-----------------%
  %---------------------------%
  
end
%-------------------------------------%

%-------------------------------------%
%-subfunctions

%---------------------------%
function cbeta = get_beta(SPM, XYZ)
%-get average of estimate betas
pwdir = pwd;
cd(SPM.swd) % you never know with relative paths of SPM

beta  = spm_get_data(SPM.Vbeta, XYZ');
beta = mean(beta, 2);
cbeta = SPM.xCon(1).c'*beta;
cbeta = mean(cbeta);

cd(pwdir)
%---------------------------%

%---------------------------%
function output = compare_con(SPM, ROI)
%-check contrasts in LC

%-------%
%-ROI in voxel space
XYZ = SPM.VM.mat \ [ROI'; ones(1,size(ROI,1))];
XYZ = XYZ(1:3,:)';
%-------%

output = '';

%-----------------%
for v = 1:size(XYZ,1) % for each voxel in ROI
  
  for b = 1:size(SPM.xX.X,2)
        
    clear allval
    
    conall = find(SPM.xX.X(:,b));
    for j = 1:numel(conall)
    confile = SPM.xY.P{conall(j)};
      
      val  = spm_get_data(confile, XYZ(v,:)');
      allval(j) = mean(val);
      
    end
    
    [~, P] = ttest(allval);
    if P < .1
    output = sprintf('%s   ROI [% 3.f % 3.f % 3.f] - base n. %1.f: % 7.3f, p-value% 3.4f\n', ...
      output, ROI(v,1), ROI(v,2), ROI(v,3), b, mean(allval), P);
    end
  end
end
%-----------------%
%---------------------------%
%-------------------------------------%