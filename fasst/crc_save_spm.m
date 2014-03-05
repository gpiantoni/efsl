function save_spm(prefix,files,d,data)

% FORMAT save_spm(prefix,files,d,data)
% Save EEG data in SPM format
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_save_spm.m 184 2010-05-17 21:10:31Z christophe $

% If, by any chance, the Meeg object is still in the structure, I will get
% rid of it before saving.

if isfield(d,'Radc')
    d = crc_spm5tospm8(files);
end

if isfield(d,'Dmeg')
    d = rmfield(d,'Dmeg');
end


[a,b] = fileparts(files(1,:));
[path name ext]=fileparts(files(1,:));

cleaned_data=d;
cleaned_data.data.fnamedat=[prefix name '.dat'];

if exist(strcat(path,filesep,cleaned_data.data.fnamedat),'file')
  delete(strcat(path,filesep,cleaned_data.data.fnamedat))
  disp(['!!! WARNING : EXISTING ' cleaned_data.data.fnamedat ' FILE OVERWRITTEN !!!' ])
end

fpd_clean = fopen(strcat(path,filesep,cleaned_data.data.fnamedat), 'a'); % 'a' append

%write the data in file .dat
fwrite(fpd_clean, data, 'float32');
cleaned_data.data.scale = ones(size(data, 1), 1);

fclose(fpd_clean);

cleaned_data.fname = [prefix name '.mat'];
cleaned_data.Nsamples = size(data,2);

D = cleaned_data;
D.data.y = file_array(fullfile(D.path, D.data.fnamedat), ...       % fname     - filename
                    [length(D.channels) D.Nsamples size(D.trials,1)],...  % dim       - dimensions (default = [0 0] )
                    spm_type('float32'), ...               % dtype     - datatype   (default = 'uint8-le')
                    0, ...                                  % offset    - offset into file (default = 0)
                    D.data.scale);                               % scl_slope - scalefactor (default = 1)
D.data.datatype = 'float32';
save(fullfile(a,D.fname),'D')
