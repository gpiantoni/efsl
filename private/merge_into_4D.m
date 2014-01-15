function single_nii_4D = merge_into_4D(list_nii_3D)
%MERGE_INTO_4D merge multiple 3D files into one 4D nii
%
% Parameters
% ----------
% list_nii_3D : cell of str
%   list with the names of the 3D files. If it ends with ',1', ',1' will be
%   appended to the output.
% 
% Returns
% -------
% single_nii_4D : str
%   name of the 4D nii file. 
% 
% see also split_into_3D

comma = strfind(list_nii_3D{1}, ',');
if ~isempty(comma)
  use_vol_idx = true;
  for i = 1:numel(list_nii_3D)
    list_nii_3D{i} = list_nii_3D{i}(1:comma-1);
  end
else
  use_vol_idx = false;
end

single_nii_4D = [list_nii_3D{1}(1:end-10) '.nii'];
spm_file_merge(list_nii_3D, single_nii_4D);
    
for i = 1:numel(list_nii_3D)
  delete(list_nii_3D{i})
end

if use_vol_idx
  single_nii_4D = [single_nii_4D ',1'];
end    
