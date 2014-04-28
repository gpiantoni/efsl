function list_nii_3D = split_into_3D(single_nii_4D)
%SPLIT_INTO_3D split one 4D nii into multiple 3D files
%
% Parameters
% ----------
% single_nii_4D : str
%   name of the 4D nii file. If it ends with ',1', ',1' will be appended to
%   the output.
% 
% Returns
% -------
% list_nii_3D : cell of str
%   list with the names of the 3D files.
% 
% see also merge_into_4D

comma = strfind(single_nii_4D, ',');
if ~isempty(comma)
  use_vol_idx = true;
  single_nii_4D = single_nii_4D(1:comma-1);
else
  use_vol_idx = false;
end

nii_3D = spm_file_split(single_nii_4D);
delete(single_nii_4D)
delete([single_nii_4D(1:end-3) 'mat'])

list_nii_3D = {};
for i = 1:numel(nii_3D)
  if use_vol_idx
    list_nii_3D{i} = [nii_3D(i).fname ',1'];
  else
    list_nii_3D{i} = nii_3D(i).fname;
  end
end