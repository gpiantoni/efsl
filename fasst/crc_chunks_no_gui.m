function crc_chunks_no_gui(opt)
%CRC_CHUNKS_NO_GUI make chunks based on crc_chunks, without gui
%
% I (GP) want to use the functions in crc_chunks of FAST, but it calls up
% the GUI all the time. Here we just copy the functions that are relevant
% to use (only using markers for example). The important thing is that the
% end result (the chunk) is identical to the output of crc_chunks.m)
% 
% OPT
%  .eeg_file : path to .mat file

handles = [];
handles.Dmeg = spm_eeg_load(opt.eeg_file);

handles.clocktime = false;

% line 525
evt=events(handles.Dmeg);

% gp
Begmarknum = find([evt.value] == opt.begmark);
Endmarknum = find([evt.value] == opt.endmark);


% line 565
Begpts = round(max(min(evt(Begmarknum).time*fsample(handles.Dmeg),nsamples(handles.Dmeg)),1));

% line 604
Endpts=round(max(min(evt(Endmarknum).time*fsample(handles.Dmeg),nsamples(handles.Dmeg)),1));

% line 614
process_chunk(handles.Dmeg,Begpts,Endpts, handles.clocktime, opt.idx_chk)

function process_chunk(Dmeg,Begpts,Endpts, clocktime, numchk_GP)
%GP: copied directly from crc_chunks, apart from numchnk_GP, which comes
%from argument (so we keep track of the name of the chunk)

d = struct(Dmeg);

try
    winsize=median([d.other.CRC.score{3,:}]);

    firstwin=floor(Begpts/(fsample(Dmeg)*winsize))+1;
    Begpts=((firstwin-1)*winsize+1/fsample(Dmeg))*fsample(Dmeg);

    lastwin=ceil(Endpts/(fsample(Dmeg)*winsize));
    Endpts=lastwin*winsize*fsample(Dmeg);

    for nsc=1:size(d.other.CRC.score,2)
        d.other.CRC.score{1,nsc}=d.other.CRC.score{1,nsc}(firstwin:lastwin);
        d.other.CRC.score{4,nsc}=[1/fsample(Dmeg) Endpts-Begpts-1/fsample(Dmeg)];

        d.other.CRC.score{5,nsc}=d.other.CRC.score{5,nsc}-Begpts;
        itemtosupr=unique(mod(find(or(d.other.CRC.score{5,nsc}<0,d.other.CRC.score{5,nsc}>Endpts-Begpts)),size(d.other.CRC.score{5,nsc},1)));
        idxkept=setdiff(1:size(d.other.CRC.score{5,nsc},1),itemtosupr);
        d.other.CRC.score{5,nsc}=d.other.CRC.score{5,nsc}(idxkept,:);
        
        d.other.CRC.score{6,nsc}=d.other.CRC.score{6,nsc}-Begpts;
        itemtosupr=unique(mod(find(or(d.other.CRC.score{6,nsc}<0,d.other.CRC.score{6,nsc}>Endpts-Begpts)),size(d.other.CRC.score{6,nsc},1)));
        idxkept=setdiff(1:size(d.other.CRC.score{6,nsc},1),itemtosupr);
        d.other.CRC.score{6,nsc}=d.other.CRC.score{6,nsc}(idxkept,:);
        
        d.other.CRC.score{7,nsc}=d.other.CRC.score{7,nsc}-Begpts;
        itemtosupr=unique(mod(find(or(d.other.CRC.score{7,nsc}<0,d.other.CRC.score{7,nsc}>Endpts-Begpts)),size(d.other.CRC.score{7,nsc},1)));
        idxkept=setdiff(1:size(d.other.CRC.score{7,nsc},1),itemtosupr);
        d.other.CRC.score{7,nsc}=d.other.CRC.score{7,nsc}(idxkept,:);
        
    end
end

Maxmem = 100*1024^2 ; % Maxmem in byte

% Assuming data are in 32 bits/ 4 byte
% Definition of the precision

% if strfind(d.data.y.dtype,'64')
%     prec = 8;
% elseif strfind(d.data.y.dtype,'32')
%     prec = 4;
% elseif strfind(d.data.y.dtype,'16')
%     prec = 2;
% elseif strfind(d.data.y.dtype,'8')
%     prec = 1;
% else
%     prec = 8;
% end

Maxelts=floor(Maxmem/(nchannels(Dmeg)*4));

Nmbchk=ceil((Endpts-Begpts)/Maxelts);

numchk = 1;

while exist([path(Dmeg) '\chk' num2str(numchk) '_' fname(Dmeg)],'file')
    numchk=numchk+1;
end

prefix = ['chk' num2str(numchk_GP) '_'];

% Save first chunk.

h = waitbar(0,'Please wait...');

cleaned_data = Dmeg(:,Begpts:min([Begpts+Maxelts nsamples(Dmeg) Endpts]));


% Update structure and save file

% handling the starting time information if available

if clocktime  
    hour = d.other.info.hour;
    date = d.other.info.date;

    newbeg = datevec(datenum([date hour]) + datenum([0 0 0 crc_time_converts(Begpts/fsample(Dmeg))]));

    d.other.info.hour = newbeg(4:6);
    d.other.info.date = newbeg(1:3);
end

for mm=1:length(d.trials.events)
    d.trials.events(mm).time = d.trials.events(mm).time - Begpts/(d.Fsample) ;
end

% Look for events beginning before the beginning of the new file and
% suppressing them. (The new file begins at the first corrected scans, i.e.
% the (1+Nscig)th scan.)
tokeep = find(~(or([d.trials.events.time]<0,[d.trials.events.time]>(Endpts-Begpts)/fsample(Dmeg))));

d.trials.events=d.trials.events(tokeep);

% Handle Chunking

file = [d.path ,filesep, d.fname];

crc_save_spm(prefix,file,d,cleaned_data);

file2 = [d.path ,filesep, prefix d.fname];


for ii=2:Nmbchk
    string = ['Please wait... ' num2str(100*((ii-1)/Nmbchk)) ' %'];
    waitbar((ii-1)/Nmbchk,h,string)
    cleaned_data = Dmeg(:,Begpts+Maxelts*(ii-1)+1:min([Begpts+Maxelts*ii nsamples(Dmeg) Endpts]));
    load(file2)
    if isfield(D,'Radc')
        D = crc_spm5tospm8(file);
    end
    d2 = D;
    crc_append_spm(file2,d2,cleaned_data)
    d2.Dmeg=meeg(D);
end

string = ['Please wait... ' num2str(100*1) ' %'];

waitbar(1,h,string)


close(h)
