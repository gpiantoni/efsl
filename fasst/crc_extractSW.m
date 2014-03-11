function [SW, cleanSW]=crc_extractSW(Ds,stages)

%function to extract SWS stages of whole nights data sets.
%
%inputs: scored whole night recording and scores to extract (vector or
%double)
%
%outputs: 
%  - SW: new file with prefix SWS_ containing only SWS stages, without
%        artefacts or arousals
%  - cleanSW: vector containing the time point position in the original 
%        file of each selected data point
%  
%__________________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre
%
% Written by J. Schrouff & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_extractSW.m 203 2010-05-26 15:08:28Z jessica $

if nargin<1
    file = spm_select;
    try
        Ds=spm_eeg_load(file);
    catch
        load(file);
        D=crc_spm5tospm8(D);
        Ds=meeg(D);
    end
end
if nargin<2
    stages=[3 4];
end

if ~isfield(Ds.CRC, 'score')
    warning('File not scored: please score sleep using FAST!')
    return
end

sc=Ds.CRC.score;
rate=Ds.fsample;

%get length of time window used to score
lw=sc{3};
sww=[];
%get windows scored as wished by the user
for i=1:size(stages,2)
    sw=find(sc{1}==stages(i));
    sww=union(sww,sw);
end

%build vector of corresponding time bins
swb=zeros(lw*rate, size(sww,2));

for i=1:size(sww,2)
    swb(:,i)=((sww(i)-1)*lw*rate:sww(i)*lw*rate-1)';
end
swb=swb(:);

%get artefacts and arousals in considered windows
%artefacts
artb=[];
if ~isempty(sc{5})
    for i=1:size(sc{5},1)
        art=sc{5}(i,1)*rate:sc{5}(i,2)*rate;
        if ~isempty(intersect(art,swb))
            artb=[artb, art];
        end
    end
end
%arousals
arou=[];
if ~isempty(sc{6})
    for i=1:size(sc{6},1)
        aro=sc{6}(i,1)*rate:sc{6}(i,2)*rate;
        if ~isempty(intersect(aro,swb))
            arou=[arou, aro];
        end
    end
end

if ~isempty(artb) && ~isempty (arou)
    disc=unique(artb, arou);
elseif ~isempty(artb) && isempty (arou)
    disc=artb;
elseif isempty(artb) && ~isempty (arou)
    disc=arou;
else
    disc=[];
end
    

%discard them from selected windows
cleanSW=setdiff(swb,disc);

clear swb artb arou disc sww

% %create new data structure
% D = clone(Ds, ['SW_' fnamedat(Ds)], [Ds.nchannels size(cleanSW,1) Ds.ntrials]);
% 

%save only clean SWS stages, avoiding 'out of memory' errors
szd=size(cleanSW,1)*Ds.nchannels*8;
memsz  = 1024^2*200;%1/5*spm('Memory');
numblocks=ceil(szd/memsz);
sizblocks=ceil(size(cleanSW,1)/numblocks);
SWD=struct(Ds);

%get first block to init new structure
stopdat=min(size(cleanSW, 1),sizblocks);
swdat=Ds(:,cleanSW(1:stopdat));
prefix='SW_';

crc_save_spm(prefix,[Ds.path, filesep,Ds.fname],SWD,swdat);
clear swdat

%append other blocks
SW=spm_eeg_load([Ds.path,filesep,prefix,Ds.fname]);
SWD=struct(SW);
idat=sizblocks+1;
for i=2:numblocks
    stopdat=min(size(cleanSW, 1),idat+sizblocks-1);
    swdat=Ds(:,cleanSW(idat:stopdat));
    if ~isempty(swdat)
        crc_append_spm(SW.fname,SWD,swdat);
    end
    idat=idat+sizblocks;
    clear swdat SW SWD
    SW=spm_eeg_load([Ds.path,filesep,prefix,Ds.fname]);
    SWD=struct(SW);
end

save(SW);







