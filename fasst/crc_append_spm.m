function crc_append_spm(files,d,data)
% FORMAT crc_append_spm(files,d,data)
%
% Append the data in matrix "data" into the *dat file of the file "files"
% and save the result in "prefixfiles.mat/dat".
%
% We are assuming that the data are in float32 format !
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_append_spm.m 184 2010-05-17 21:10:31Z christophe $

cleaned_data = d;

[a,b] = fileparts(files(1,:));

fpd_clean = fopen(fullfile(a, cleaned_data.data.fnamedat), 'a'); % 'a' append

% Write the data in file .dat
fseek(fpd_clean,0,'eof');
fwrite(fpd_clean, data, 'float32');
cleaned_data.data.scales = ones(size(data, 3-2), 1);
    
fclose(fpd_clean);

cleaned_data.Nsamples = d.Nsamples+size(data,2);

D = cleaned_data;
D.data.y = file_array( ...
    fullfile(D.path, D.data.fnamedat), ... % fname - filename
    [length(D.channels) D.Nsamples size(D.trials,1)],...  % dim - dimensions (default = [0 0] )
    spm_type('float32'), ... % dtype - datatype   (default = 'uint8-le')
    0, ... % offset - offset into file (default = 0)
    D.data.scale); % scl_slope - scalefactor (default = 1)
D.data.datatype = 'float32';
save(fullfile(a,D.fname),'D')
