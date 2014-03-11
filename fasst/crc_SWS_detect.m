function [D]= crc_SWS_detect(handles)

%detects Slow Waves Sleep (and delta waves)on EEG data. The criteria are
%based on Massimini's but were softened to detect all waves. The data are 
%first loaded and stages 3 and 4 are extracted if the file is scored using 
%FAST. If the file is not scored, we supposed that it contains only stages
%3 or 4. The data are then filtered (0.2-4Hz) if needed. The detection is
%achieved on 4 Regions Of Interest (frontal, left central, right central,
%centro-parietal, which can be either automatically created or manually) 
%and trajectory of each wave is then detected on the whole scalp. 
%
% The criteria used are magnitude criteria in microV ( minimum negative peak amplitude
% -40 for delta and -80 for SWS, minimum total magnitude: 75 for delta and
% 140 for SWS) and duration criteria in ms (duration of negative peak:
% between 250 and 1250 ms, time between up zero crossing and positive 
%peak : maximum 2000 ms). Another criterion was added on the slope between
%the negative and the positive peaks (criteria of  minimum percentile 90).
%
% The outputs are new structures in D(data).other:
% D.other.SW where SW is a structure with:
%SW=struct('down', [], ...  %index of downward zero crossing
%     'up', [], ...          % index of upward zero crossing
%     'negmax', [], ...      % index of maximum negativity using Massimini criteria
%     'negmax_tp', [], ...   % negmax in time points
%     'posmax_tp', [], ...   % index of maximum positivity
%     'posmax', [], ...      % posmax in time points
%     'upstart', [], ...     % start of the upstate
%     'upend', [], ...       % end of upstate
%     'maxslope', [], ...    % maximum of slope index in the upswing
%     'channels', [], ...    % for each wave, channel numbers respecting the criteria in E.data
%     'electrodes', [], ...  % for each wave, channel names respecting the criteria in E.data
%     'delays', [], ...      % for each wave, delay of the minimum in the SW, for each channel
%     'uponset', [], ...     % for each wave, the onset of the up state in terms of scan
%     'amplitude', [],...    % peak to peak magnitude
%     'code',[],...          % code of the wave
%     'DATA4ROI',[],...      % mean of the data on 4 regions of interest
%     'neg_slope',[],...     % maximum of slope between DZC and negmax
%     'negmax_TR',[]);       % negmax in TR, computed after the beginning of the scan
%
% If the user chose to review the detected waves, another structure will be
% created containing only 'accepted' waves: D.other.goodSW
%
%The code of the wave is 3'(number of the first electrode detecting it)' for
%SWS and 4'(number of the first electrode detecting it)' for delta waves.
%All fields are in seconds (or in time points when said explicitly).
%
%The field .events.origin_count counts the number of waves 'starting' 
% with each electrode and the total number of waves detected by each 
% electrode.
%
% The field D.trials.events was also modified to contain the type 'delta'
% for a delta wave and 'SW' for a slow wave, the value equals to the code of 
%the wave, the time corresponding to the maximum power of the negative peak
%(in seconds) and the duration corresponding to the maximum delay.
%
%The possibility to review all waves is given and allows the user to
%reject some waves. The 'accepted' waves will be saved in D.other.goodSW,
%separately from the automatically detected waves.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%This program was elaborated with the help of Pierre Maquet and Christophe
%Phillips who have kindly shared routines and parts of code.
%Jessica Schrouff, 23/08/2007, CRC Ulg
% last modified 19/01/2009 for SPM8
% last modified 22/10/2009 for FAST toolbox
% last modified 11/05/2010 for chunked files and default electrodes
% positions
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by J. Schrouff & C. Phillips, 2009.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_SWS_detect.m 218 2010-05-31 14:36:43Z jessica $

global gindex
close all

% loading files
%--------------------------------------------------------------------------
if ~nargin
    file = spm_select(1, 'mat', 'Select cleaned EEG file','',pwd,'.*');
    [a,b] = fileparts(file);
    Dss = spm_eeg_load(strcat(a,filesep,b,'.mat'));
else
    Dss = spm_eeg_load(handles.fname);
end
close all
XV=struct(Dss);

def=crc_get_defaults();  % GP
def.swsd.param1 = def.swsd;  % GP
csw=[];
start=0;
%extract periods of interest
if handles.analyse==3
    if isfield(Dss.CRC, 'score') && size(XV.other.CRC.score,1)>4
        if ~isempty(handles.stagesw)
            [Ds, cleanSW]=crc_extractSW(Dss, handles.stagesw);
        else
            [Ds, cleanSW]=crc_extractSW(Dss, def.swsd.stagesw);
        end
        csw=cleanSW/fsample(Dss);
    end
elseif handles.analyse==2
    Ds=Dss;
    disp('Warning: assuming file contains only SWS stages')
else
    crc_process_chunk(Dss,handles.Begpts,handles.Endpts,1);
    numchk = 1;
    while exist([path(Dss) '\chk' num2str(numchk) '_' fname(Dss)],'file')
        numchk=numchk+1;
    end
    prefix = ['chk' num2str(numchk) '_'];
    Ds=spm_eeg_load([path(Dss),filesep,prefix,fname(Dss)]);
    start=handles.Begpts/fsample(Ds); %relative start in s
end
D=struct(Ds);

%consider joint fMRI-EEG case
startscan=0;
%get time of scan start w.r.t EEG start 
if handles.fmri && ~isempty(handles.TR) && ~isempty(handles.marker)
    for i=1:size(D.trials.events,2)
        if strcmpi(D.trials.events(i).value,handles.marker)
            startscan=D.trials.events(i).time;
            break
        end
    end
    if startscan==0
        disp('Marker not found in EEG events structure!!')
        cont=spm_input('Continue?',1,'y|n',[1 0],0);
        close gcf
        if cont==0
            return
        end
    end
    %takes it into account for timing of waves
    if ~isempty(csw)
        csw=csw-startscan;
    else
        start=start-startscan;       
    end
end

param=struct('t',0:(1/D.Fsample):(size(D.data.y,2)/D.Fsample-1/D.Fsample),...
            'SWmAmpl',[-40 -80 75 140]);

clc
D_channels_eeg=find(strcmpi('eeg',{D.channels.type}));

%Checking montage
%------------------------------------------------------------------
%if no re-reference field in D.history, then warning and ask to continue or
%not

disp('---- Checking reference----')
if handles.reref==0
    disp('WARNING: no re-referencing found!-- please re-reference using spm8 (spm_eeg_montage)')
    cont=spm_input('Continue?',1,'y/n',[1,0],1);
else
    disp('re-referencing achieved....ok')
    cont=1;
end
if cont~=1
    return
end

%filtering data
%--------------------------------------------------------------------------
%Use of a lowpass filter (cutoff rate asked in gui) followed by a highpass
%filter. Filtering is achieved during detection to avoid writing new data 
%on the hard drive.
%For Slow Waves Sleep detection, we recommand lowpass filtering at 4Hz and 
%highpass filtering at 0.25.

if ~isempty(handles.highfc)&& ~isempty(handles.lowfc)
    fqcut=[handles.highfc, handles.lowfc];
else
    fqcut=[def.swsd.highfc, def.swsd.lowfc];
end
disp('-----filtering data-----')
args=[];
dd=0;
for i=1:size(D.history,2)
    if strcmpi('spm_eeg_filter',D.history(i).fun)
        dd=dd+1;
        args=[args, D.history(i).args.filter];
    end
end
vv=0;        
PHz=zeros(1,size(args,2));
if ~isempty(args)
    for i=1:size(args,2)
        PHz(i)=args(i).PHz(1);
    end
    fc=sort(PHz);
    if fc(1)<0.3 || fc(2)<10 || fc(2)>3
        vv=1;
    end
end

fhc = fqcut(1)/(D.Fsample/2);
flc = fqcut(2)/(D.Fsample/2);
order=def.swsd.butterorder;
[b1,a1]=butter(order,fhc,'high');
[b2,a2]=butter(order,flc,'low');

%creating non overlapping regions of interest
%--------------------------------------------------------------------------


disp('---creating non overlapping ROIs---') 
DATA4ROI = struct('data',[],'channels',[],'rate',[]);

auto=handles.roisel;
h = waitbar(0,'Please wait...');
if auto==1
    %automatic selection of ROI based on electrodes position in
    %crc_electrodes.mat (2D)
    load CRC_electrodes.mat
    el_names_temp=cell(size(names,2),1);
    for i=1:size(names,2)
        el_names_temp{i}=upper([names{i}]);
    end
    [origin_count_dat] = pm_origin_count(D);
    [dumb1,dumb2,index2]=intersect(origin_count_dat(:,1),el_names_temp);
    eeg_chan=index2(find(crc_types(index2)>-2));
    ROI_centers=def.swsd.ROI_centers;
    for i=1:size(ROI_centers,1)
        pos_eeg_chan=pos(eeg_chan,:);
        pos_cent=ones(size(pos_eeg_chan,1),2);
        pos_cent(:,1)=ROI_centers(i,1);
        pos_cent(:,2)=ROI_centers(i,2);
        dist_pos=((pos_cent(:,1)-pos_eeg_chan(:,1)).^2+(pos_cent(:,2)-pos_eeg_chan(:,2)).^2).^0.5;
        roi_chan=find(dist_pos<=0.1);
        [fin_roi,ind1]=intersect(upper(origin_count_dat(:,1)),upper(names(eeg_chan(roi_chan))));
        tmp = zeros(1,size(D.data.y,2));
        div=1;
        
        for j=1:length(roi_chan)
            importdat=D.data.y(ind1(j),:);
            if vv~=1
                importdat = filtfilt(b1, a1, importdat);
                importdat = filtfilt(b2, a2, importdat);
            end
            tmp(1,:) = tmp(1,:)*(1-(1/div))+importdat*(1/div);
            div=div+1;
            
        end
        DATA4ROI.data= [DATA4ROI.data;tmp];
        string = ['Please wait... ' num2str(100*(i/size(ROI_centers,1))) ' %'];
        waitbar((i)/size(ROI_centers,1),h,string)
    end
elseif auto==0
    numroi=handles.numroi;
    name_roi=handles.name_roi;
    for iroi = 1:numroi %regions of interest
        tmp = zeros(1,size(D.data.y,2));
        div=1;
        sel_ROI=handles.sel_ROI;
        for iselroi=1:size(sel_ROI{iroi},2)
            importdat = D.data.y(sel_ROI{iroi}(iselroi),:);
            if vv~=1
                importdat = filtfilt(b1, a1, importdat);
                importdat = filtfilt(b2, a2, importdat);
            end
            tmp(1,:) = tmp(1,:)*(1-(1/div))+importdat*(1/div);
            div=div+1;
            
        end
        DATA4ROI.data= [DATA4ROI.data;tmp];
        string = ['Please wait... ' num2str(100*(iroi/numroi)) ' %'];
        waitbar((iroi)/numroi,h,string)
    end
end
string = ['Please wait... ' num2str(100*1) ' %'];
waitbar(1,h,string)
close(h)

% Completing DATA4ROI structure
DATA4ROI.roisel=auto;
if auto==0
    DATA4ROI.nameroi=name_roi;
else
    DATA4ROI.nameroi=[];
end
DATA4ROI.channels = D.channels;
DATA4ROI.rate = D.Fsample;
clear  tmp ichannel ielec iroi 

               
%detecting criteria over the ROIs
%--------------------------------------------------------------------------

disp(['---detecting criteria over the ROIs---'])

%structure for each wave

SW=struct('down', [], ...  %index of downward zero crossing
    'up', [], ...          % index of upward zero crossing
    'negmax', [], ...      % index of maximum negativity using Massimini criteria
    'negmax_tp', [], ...   % negmax in time points
    'posmax', [], ...      % index of maximum positivity
    'posmax_tp', [], ...   % posmax in time points
    'upstart', [], ...     % start of the upstate
    'upend', [], ...       % end of upstate
    'maxslope', [], ...    % maximum of slope index in the upswing
    'channels', [], ...    % for each wave, channel numbers respecting the criteria in E.data
    'electrodes', [], ...  % for each wave, channel names respecting the criteria in E.data
    'delays', [], ...      % for each wave, delay of the minimum in the SW, for each channel
    'uponset', [], ...     % for each wave, the onset of the up state in terms of scan
    'amplitude', [],...    % peak to peak magnitude
    'neg_slope', [],...     % maximum of negative slope
    'code',[],...          % code of the wave: 4 for delta and 3 for SWS followed by the number of the channel
    'negmax_TR',[]);       % negmax in TR, computed after the beginning of the scan

param1=def.swsd.param1;
countwaves=0;
            
for idataf = 1:size(DATA4ROI.data,1) % looping over the 4 ROIs
    disp(['ROI number ' num2str(idataf) ])
    F = DATA4ROI.data(idataf,:);
    [SW,countwaves,D] = find_SW(F,DATA4ROI,D,param1,countwaves,SW);
    disp(['Total number of waves detected:' num2str(countwaves)])
end

if ~(size(SW,2)>1)
    disp('-- No Slow Waves detected in this data set --')
    return
end

%detecting waves on each channel - establishing travel directory
%--------------------------------------------------------------------------
%detects waves on all electrodes and establishes the trajectory followed by
%the wave on the scalp.

disp(['----establishing traveling direction----'])

ONSETS= struct('channel',[],'onsets',[]);
SWS_count=0;
origin_count=pm_origin_count(D);
delta_count=0;
for jchannel = 1:size(origin_count,1)
    ONSETS(jchannel).channel = origin_count{jchannel,1};
end
evtsw=struct('type',[],'value',[],'duration',[],'time',[],'offset',[]);

%Define additional indexes to load in order to filter properly and
%therefore avoid windowing effects

[h1,t1]=impz(b1,a1);
[h2,t2]=impz(b2,a2);
addpts=max(size(t2,2),size(t1,1));
addpts=min(addpts,20000);
h = waitbar(0,'Please wait...');

for iwav = 1:size(SW,2)                          % for each detected wave

    %     Mind that the index cannot be <0
    sampleindex =round(SW(iwav).negmax_tp-param1.SWlength(4)*D.Fsample/1000:SW(iwav).negmax_tp+param1.SWlength(4)*D.Fsample/1000);
    extended_idx=[sampleindex(1)-addpts:sampleindex(1)-1, sampleindex, sampleindex(end)+1:sampleindex(end)+addpts];
    markervect=[zeros(1,addpts) ones(1,length(sampleindex)) zeros(1,addpts)];
    
    tokeep = find(and(extended_idx > 0, extended_idx < D.Nsamples));

    extended_idx=extended_idx(tokeep);
    markervect=markervect(tokeep);
    
    sample = D.data.y(D_channels_eeg,extended_idx);             % display 10 s
    if vv~=1 %filter data if necessary
        for zchan=1:size(sample,1)
            sample(zchan,:) = filtfilt(b1, a1, sample(zchan,:));
            sample(zchan,:) = filtfilt(b2, a2, sample(zchan,:));
        end
    end

    recupidx=find(markervect==1);
    sample=sample(:,recupidx);
    
    
    rds=sample(:,round(size(sample,2)/2)-20:round(size(sample,2)/2)+20);
    power=sum(rds.^2);  %time = max of power
    [ipower,jpower]=max(power);
    SW(iwav).negmax_tp=SW(iwav).negmax_tp+(jpower-round(size(rds,2)/2));
    SW(iwav).start=SW(iwav).negmax_tp;
    SW(iwav).negmax=SW(iwav).negmax_tp*1000/D.Fsample;
    [sampmin,sampos] = min(sample');
    [sortedsamp, samporder] = sort(sampos);          % sortedsamp = ordered delays; samporder = ranking order

    %     delete shallow and 'up' channels
    deriv = diff(sample(samporder,:)')';
    if SW(iwav).code == 100
        samp_delta =(sampmin(samporder) < param.SWmAmpl(1));
        sampindex = samp_delta .* (mean(deriv(:,1:round(size(sample,2)/2))') < 0);
        delta_count=delta_count+1;% modified for delta waves
    elseif SW(iwav).code == 101
        samp_SWS = (sampmin(samporder) < param.SWmAmpl(2));
        sampindex = samp_SWS .* (mean(deriv(:,1:round(size(sample,2)/2))') < 0);
        SWS_count=SWS_count+1;
    end
    all_names={D.channels(:).label};
    channels_eeg = all_names(D_channels_eeg);
    sortdelays = sortedsamp(find(sampindex));
    sortchannel = samporder(find(sampindex));
    SW(iwav).channels = sortchannel; %channels respecting magnitude and slope criteria, ordered by delays
    SW(iwav).electrodes = channels_eeg(sortchannel); % corresponding electrodes
    SW(iwav).delays = sortdelays; % increasing delays (of minimum) 
    
    %define new codes
    if SW(iwav).code == 101
        if length(num2str(SW(iwav).channels(1)))==1
            a=['30',num2str(SW(iwav).channels(1))];
        elseif length(num2str(SW(iwav).channels(1)))==2
            a=['3',num2str(SW(iwav).channels(1))];
        end
        evtsw(iwav).type='SW';
    elseif SW(iwav).code == 100
        if length(num2str(SW(iwav).channels(1)))==1
            a=['40',num2str(SW(iwav).channels(1))];
        elseif length(num2str(SW(iwav).channels(1)))>=2
            a=['4',num2str(SW(iwav).channels(1))];
        end
        evtsw(iwav).type='delta';
    end
    evtsw(iwav).value=str2double(a);
    SW(iwav).code=str2double(a);
    evtsw(iwav).time=SW(iwav).negmax_tp/D.Fsample;
    evtsw(iwav).duration=max(sortdelays); 
   
    % count the SW channel start
    
    for kelec = 1:size(origin_count,1)
        if strcmpi((origin_count{kelec,1}),char(SW(iwav).electrodes(1)))
            origin_count{kelec,2} = origin_count{kelec,2}+1;
        end
        for i=1:size(SW(iwav).electrodes,2)
            if strcmpi((origin_count{kelec,1}),char(SW(iwav).electrodes(i)))
                origin_count{kelec,3} = origin_count{kelec,3}+1;
            end
        end
    end
    
    %get time in TR if joint EEG-fMRI acquisition
    if handles.fmri && startscan~=0
        SW(iwav).negmax_TR=(SW(iwav).negmax+start)/handles.TR;
    else
        SW(iwav).negmax_TR=NaN;
    end
    string = ['Please wait... ' num2str(100*(iwav/size(SW,2))) ' %'];
    waitbar((iwav)/size(SW,2),h,string)
    
end
string = ['Please wait... ' num2str(100*1) ' %'];
waitbar(1,h,string)


%save file with data of interest
D.other.SW = SW;
D.other.origin_count=origin_count;
D.other.DATA4ROI = DATA4ROI;
D.trials.events=[D.trials.events, evtsw];

%compute the timing of waves w.r.t. original data file
for iwav = 1:size(SW,2)
    if handles.analyse==3
        start=csw(SW(iwav).start);
        SW(iwav).down = csw(SW(iwav).down/1000*D.Fsample)*1000;         %position of downward zero crossing
        SW(iwav).up = csw(SW(iwav).up/1000*D.Fsample)*1000;             %position of upward zero crossing
        SW(iwav).negmax = csw(SW(iwav).negmax_tp)*1000;                 %negative peak position in ms
        SW(iwav).negmax_tp= csw(SW(iwav).negmax_tp);                    %negative peak position in time points
        SW(iwav).posmax = csw(SW(iwav).posmax_tp)*1000;                 %positive peak position in ms
        SW(iwav).posmax_tp= csw(SW(iwav).posmax_tp);                    %positive peak position in time points
        SW(iwav).upstart = csw(SW(iwav).upstart/1000*D.Fsample)*1000;   %begin of upstate
        SW(iwav).upend = csw(SW(iwav).upend/1000*D.Fsample)*1000;       %end of upstate
        %get time in TR if joint EEG-fMRI acquisition
        if handles.fmri && startscan~=0
            SW(iwav).negmax_TR=(SW(iwav).negmax)/handles.TR;
        else
            SW(iwav).negmax_TR=NaN;
        end
    else
         SW(iwav).down = SW(iwav).down+start*1000;                                
         SW(iwav).up = SW(iwav).up+ start*1000;                                   
         SW(iwav).negmax = SW(iwav).negmax+ start*1000;               
         SW(iwav).negmax_tp= SW(iwav).negmax_tp+start;               
         SW(iwav).posmax = SW(iwav).posmax+ start*1000;  
         SW(iwav).posmax_tp= SW(iwav).posmax_tp+start;                
         SW(iwav).upstart = SW(iwav).upstart+ start*1000;  
         SW(iwav).upend = SW(iwav).upend + start*1000; 
         %get time in TR if joint EEG-fMRI acquisition
        if handles.fmri && startscan~=0
            SW(iwav).negmax_TR=(SW(iwav).negmax)/handles.TR;
        else
            SW(iwav).negmax_TR=NaN;
        end
    end
end
        
%save both files
XV.other.SW = SW;
XV.other.origin_count=origin_count;
XV.other.DATA4ROI = DATA4ROI;
XV.trials.events=[XV.trials.events, evtsw];

save(fullfile(D.path,D.fname),'D');
D=XV;
save(fullfile(D.path,D.fname),'D');
delete([D.path,filesep,'f',D.data.fnamedat]);
delete([D.path,filesep,'ff',D.data.fnamedat]);
close(h)

% visual check of all SWS
%--------------------------------------------------------------------------
%displays data on each region of interest near each SW and underline
%supposed SW. Subplot the scalp map where delays are interpolated.

sss = zeros(2,size(F,2));
TMP = [];
for isw = 1:size(SW,2)
    TMP = [TMP;SW(isw).negmax SW(isw).posmax];
end
TMP = round(TMP/1000*D.Fsample);
sss(1,TMP(:,1)) = def.swsd.dispscale;
sss(1,TMP(:,1)+1) = -def.swsd.dispscale;
sss(2,TMP(:,2)) = def.swsd.dispscale;
sss(2,TMP(:,2)+1) = -def.swsd.dispscale;

vc=handles.review;
% Display loop
if vc==1
    %add fileio directory of spm8 to have access to read_sens
    dspm=spm('dir');
    addpath([dspm,filesep,'external',filesep,'fileio'])
    elpos=handles.sensauto;
    if elpos==0
        if ~isempty(handles.sensfname)
            zebris_name=handles.sensfname;
        else
            zebris_name=spm_select(1, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
        end
    else
        zebris_name='CRC_electrodes.mat';
    end
    cas=handles.maps;
    accepted=[];
    for jwav = 1:size(SW,2)
        crc_SWS_mapping(D,jwav,cas,zebris_name,elpos)
        wav_ok=spm_input('Accept ?',1,'b','y|n',[1;0],1);
        if wav_ok==1
            accepted=[accepted, SW(jwav)];
        end
        D.other.goodSW=accepted;
        conti=spm_input('Show next',+1,'b','y|n',[1;0],1);
        if conti==0
            close gcf
            break
        end
    end
end

%save structure 
save(fullfile(D.path,D.fname),'D');

%--------------------------------------------------------------------------
%---- SUBFUNCTION TO DETECT WAVES CORRESPONDING TO MASSIMINI CRITERIA  ----
%--------------------------------------------------------------------------


function [SW,countwaves,E] = find_SW(F,DATA4ROI,E,param1,countwaves,SW)


% find zero crossings


F(2,:)=sign(F(1,:)); %gives the sign of data
F(3,:)=[0 diff(F(1,:))]; % gives the differential of data
F(4,:) = sign(F(3,:)); %gives the sign of differential

DZC=find(diff(F(2,:)) ==-2);
UZC=find(diff(F(2,:)) ==2);

%criterium of maximum slope index percentile 90

MSI=zeros(1,size(F,2));
MSI_plot=MSI;
MSI_plot(find(F(3,:) > prctile (F(3,:),90)))=100; % giopia
MSI=find(MSI_plot==100); 

for imsi=1:size(MSI,2)-1
    
    %find nearest MSI and DZC
    indiceDZC = find((DZC-MSI(imsi))<0);
    if ~isempty(indiceDZC)
        indiceDZC=indiceDZC(end);
        iDZC=DZC(indiceDZC);
    else
        iDZC=1;
    end
    
    
    
    %find posmin and valmin
    [valmin,indposmin]=min(F(1,iDZC:MSI(imsi)));
    posmin=iDZC+indposmin;
    negslope=min(F(3,iDZC:posmin));
    
    %find iUZC between iDZC and iDZC+SW_length(2)
    upperbound=size(F,2)-(iDZC+DATA4ROI.rate*param1.SWlength(2)/1000);
    if upperbound >0
        iUZC = find(diff(F(2,iDZC:iDZC+round(DATA4ROI.rate*param1.SWlength(2)/1000))) == 2)+iDZC; 
    else
        iUZC= find(diff(F(2,:)) ==2)+iDZC;
    end
    
    
    if ~isempty(iUZC)&& ~isempty(indiceDZC)
        iUZC=iUZC(1);
    
        
        %verification of Massimini criteria on length and magnitude of SWS
        %and delta waves
        
        %criterion on negative peak magnitude
        if ((iUZC-iDZC) <= (param1.SWlength(2)*DATA4ROI.rate/1000) && (param1.SWlength(1)*DATA4ROI.rate/1000) <=(iUZC-iDZC) )
            
            %negative peak magnitude
            if valmin <= param1.SWmAmpl(1)
                upperbound= size(F,2)-(iUZC+(DATA4ROI.rate*param1.SWlength(3)/1000));
                
                if upperbound >0
                    posmax = iUZC + find(diff(F(4,iUZC + find(F(1,iUZC:iUZC+round(DATA4ROI.rate*param1.SWlength(3)/1000)) > 0))) == -2); 
                else
                    posmax= iUZC +find(diff(F(4,iUZC+find(F(1,iUZC:size(F,2)-1 >0)))) ==-2);
                end
                
                if ~isempty(posmax)
                    start2end_up=[];
                    posmax=posmax(1);
                    valmax=F(1,posmax);
                    
                    
                 %criterion on peak to peak magnitude using mimimal
                 %criteria of -40 and 75 microV
                     if ((indiceDZC+1)< size(DZC,2)) && ((iUZC+(DATA4ROI.rate*param1.SWlength(3)/1000)) < size(F,2))

                         if DZC(indiceDZC+1)-iUZC < DATA4ROI.rate*param1.SWlength(3)/1000
                                start2end_up = iUZC + find(F(1,iUZC:MSI(imsi+1))>= 0.8*valmax); 
                         else
                                start2end_up = iUZC + find(F(1,iUZC:iUZC+round(DATA4ROI.rate*param1.SWlength(3)/1000))>= 0.8*valmax);
                         end
                     end

                     if ~isempty(start2end_up)

                         if (abs (valmax)+ abs (valmin)) >= param1.SWmAmpl(3)

                         %avoid the doubloons
                         
                          if  isempty(SW(end).negmax) || ...% nothing in SW.negmax (1st pass)
                                  ((all(abs(posmin - squeeze(cat(1,SW(:).negmax)./1000*E.Fsample))  > DATA4ROI.rate*param1.SWlength(3)/5/1000))&&...   % need 400ms between SW negativity
                                     (all(abs(posmax - squeeze(cat(1,SW(:).posmax)./1000*E.Fsample))  > DATA4ROI.rate*param1.SWlength(3)/5/1000)))     
                                
                                 %fill the SW structure being careful of
                                %delta waves (different code)
                                
                                countwaves =  countwaves+1;
                                                                
                                ms=max(F(3,posmin: posmax));
                               
                                if (valmin <= param1.SWmAmpl(2) && (abs (valmax)+ abs (valmin)) >= param1.SWmAmpl(4))
                                    code =101;
                                else 
                                    code =100;
                                end
                               
                                disp([ num2str(imsi)  ' detected possible SW - ' ...
                                    num2str(countwaves) ' SW kept in total;  at ' ...
                                    num2str(posmin) ' : ' ...
                                    num2str(valmin) ' microV; at '   ...
                                    num2str(posmax) ' : '  num2str(valmax) ' microV; '])                               
                                    
                                SW(countwaves).down = iDZC*1000/E.Fsample;                  %position of downward zero crossing
                                SW(countwaves).up = iUZC*1000/E.Fsample;                    %position of upward zero crossing
                                SW(countwaves).negmax = posmin*1000/E.Fsample;              %negative peak position in ms
                                SW(countwaves).negmax_tp= posmin;                           %negative peak position in time points
                                SW(countwaves).posmax = posmax*1000/E.Fsample;              %positive peak position in ms
                                SW(countwaves).posmax_tp= posmax;                           %positive peak position in time points
                                SW(countwaves).upstart = start2end_up(1)*1000/E.Fsample;    %begin of upstate
                                SW(countwaves).upend = start2end_up(end)*1000/E.Fsample;    %end of upstate
                                SW(countwaves).amplitude = (abs(valmax)+abs(valmin));       %total magnitude
                                SW(countwaves).uponset = SW(countwaves).negmax ;
                                SW(countwaves).maxslope = ms;
                                SW(countwaves).code= code;
                                SW(countwaves).neg_slope = negslope;
                                
                               
                            end
                        end
                    end
                end
            end
        end
    end
end

N_sw = size(SW,2);
if N_sw>1
    countwaves = N_sw;
else % check if at least one SW was found in 1st ROI(s)
    if ~isempty(SW.down)
        countwaves = N_sw;
    end
end
% countwaves= size(SW,2);

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INITIALIZE COUNT ON CHANNELS   ---------------
%--------------------------------------------------------------------------
function [origin_count] = pm_origin_count(data)
data_channels_eeg=find(strcmpi('eeg',{data.channels(:).type}));
all_names={data.channels(:).label};
origin_count=all_names(data_channels_eeg)';
for i=1:size(data_channels_eeg,2)
    origin_count{i}=upper(deblank(origin_count{i}));
    origin_count(i,2)={0};
    origin_count(i,3)={0};
end
origin_count(:,2)={0};
origin_count(:,3)={0};


function y = prctile(x,p,dim)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.


if ~isvector(p) || numel(p) == 0
    error('stats:prctile:BadPercents', ...
          'P must be a scalar or a non-empty vector.');
elseif any(p < 0 | p > 100) || ~isreal(p)
    error('stats:prctile:BadPercents', ...
          'P must take real values between 0 and 100');
end

% Figure out which dimension prctile will work along.
sz = size(x);
if nargin < 3 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),class(x));
    else
        szout = sz; szout(dim) = numel(p);
        y = nan(szout,class(x));
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    nrows = sz(dim);
    ncols = prod(sz) ./ nrows;
    x = reshape(x, nrows, ncols);

    x = sort(x,1);
    nonnans = ~isnan(x);

    % If there are no NaNs, do all cols at once.
    if all(nonnans(:))
        n = sz(dim);
        if isequal(p,50) % make the median fast
            if rem(n,2) % n is odd
                y = x((n+1)/2,:);
            else        % n is even
                y = (x(n/2,:) + x(n/2+1,:))/2;
            end
        else
            q = [0 100*(0.5:(n-0.5))./n 100]';
            xx = [x(1,:); x(1:n,:); x(n,:)];
            y = zeros(numel(p), ncols, class(x));
            y(:,:) = interp1q(q,xx,p(:));
        end

    % If there are NaNs, work on each column separately.
    else
        % Get percentiles of the non-NaN values in each column.
        y = nan(numel(p), ncols, class(x));
        for j = 1:ncols
            nj = find(nonnans(:,j),1,'last');
            if nj > 0
                if isequal(p,50) % make the median fast
                    if rem(nj,2) % nj is odd
                        y(:,j) = x((nj+1)/2,j);
                    else         % nj is even
                        y(:,j) = (x(nj/2,j) + x(nj/2+1,j))/2;
                    end
                else
                    q = [0 100*(0.5:(nj-0.5))./nj 100]';
                    xx = [x(1,j); x(1:nj,j); x(nj,j)];
                    y(:,j) = interp1q(q,xx,p(:));
                end
            end
        end
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz; szout(dim) = numel(p);
    y = reshape(y,szout);
end
% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end

                