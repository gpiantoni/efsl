function crc_SWS_mapping (D,value,choice,zebr,dim)

% FUNCTION crc_SWS_mapping
%
% maps slow waves of deep human sleep thanks to two views : 
%
%  - electrical signal averaged on 4 ROIs plotted against time in the 
%[-500ms, +1500ms] interval around the concerned (value) SW.
%
%  - 2D flattened map of the scalp showing either a hot colormap of delays
%  (choice 1) or a jet movie showing potentials on the scalp in the same
%  time interval (choice 2).
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by J. Schrouff & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_SWS_mapping.m 208 2010-05-27 13:20:02Z jessica $

if nargin<1
     file = spm_select(1, 'mat', 'Select cleaned EEG file','',pwd,'.*');
    [a,b] = fileparts(file);
    Ds = spm_eeg_load(strcat(a,filesep,b,'.mat'));
    D=struct(Ds);
 end
if nargin<2, value = spm_input ('Number of the wave to show:',1,'r',1,[]); end
close gcf
if nargin<3, choice = spm_input('Map of:',+1,'b','delays|potentials',[1;2],1); end 
close gcf
if nargin<5
    dim=spm_input('Electrodes positions:',+1,'b','file|auto',[0,1],0);
end
if nargin<4 && dim==0 
    zebr=spm_select(Inf, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
end
close gcf   
   
% close all
SW=D.other.SW;
if dim==0
    el_set_sphuni=crc_elecpos_sph(D,zebr);
else
    el_set_sphuni=[];
end

def=crc_defaults('swsd');
close gcf


%plot signal averaged on the 4 ROIs, highlighting the wave of interest.
%----------------------------------------------------------------------

%get signal of interest
for jwav = 1:size(SW,2)
    sampleindex =SW(value).start-5*D.Fsample:SW(value).start+5*D.Fsample;
    sampleindex = sampleindex(find(sampleindex < size(D.data.y,2)));
    sampleindex = round(sampleindex);
    if ~all(sampleindex > 0)
        sampleindex = sampleindex(find(sampleindex > 0));
        dispdata = D.other.DATA4ROI.data(:,sampleindex); % display 10 s
        negmax=int32(SW(value).negmax/1000*D.Fsample);
    else
        dispdata = D.other.DATA4ROI.data(:,sampleindex); % display 10 s
        negmax=int32(5*D.Fsample);
    end     
end


sss = zeros(2,size(D.data.y,2));
TMP = [];
for isw = 1:size(D.other.SW,2)
    TMP = [TMP;D.other.SW(isw).negmax D.other.SW(isw).posmax];
end

TMP = round(TMP/1000*D.Fsample);
sss(1,TMP(:,1)) = def.swsd.dispscale;sss(1,TMP(:,1)+1) = -def.swsd.dispscale;
sss(2,TMP(:,2)) = def.swsd.dispscale;sss(2,TMP(:,2)+1) = -def.swsd.dispscale;

t=0:1/D.Fsample:size(D.data.y,2)/D.Fsample-1/D.Fsample;

origin_count= pm_origin_count(D);
roi_sel=D.other.DATA4ROI.roisel;
name_roi=D.other.DATA4ROI.nameroi;
%display signal on 4 ROIs
fig=pm_display(dispdata,[10 10 800 800],75,t(sampleindex),sss(:,sampleindex), negmax,5,[-def.swsd.dispscale def.swsd.dispscale],0,[],[],roi_sel,name_roi,def); 
subplot(2,1,2)
pm_map(D.other.SW,D,origin_count,choice,value,fig,dim, el_set_sphuni,def)

%--------------------------------------------------------------------------
%---------  SUBFUNCTION TO DISPLAY DATA AND POINT WAVES  ------------------
%--------------------------------------------------------------------------

function  fig=pm_display(dispdata,position,displstep,t,sss,negmax,windowl,yLim,uiscroll,disptext,delays,roi_selec,namroi,def)

if ~nargin
    position = [10 10 800 800];
    displstep = 75;
    t = 1:size(dispdata,2);
    windowl = 10*D.Fsample;
end

%figure
fig=spm_figure('Create','1','Display of detected Slow Waves','on');
subplot(2,1,1)
axis off
hold;
WS = spm('WinScale');
position=position.*WS;
 
for jj = 1:size(dispdata,1)
    if ~isempty(disptext)
        d=(position(4))/(position(4)*10);
        axes('Position',[0.1,0.85-(jj-1)*d,0.8,d]); 
        plot(t,dispdata(jj,:)) 
        plot(delays(jj),min(dispdata(jj,:)),'o')
        axis ([t(1) t(end) yLim(1) yLim(end)]);
        str = disptext(jj);
        text(t(end)+.01,-(jj-1)*displstep,str)
        grid on
    else
        d=(position(4))/(position(4)*10);
        axes('Position',[0.1,0.85-(jj-1)*d,0.8,d]);
        plot(t,dispdata(jj,:),'LineWidth',2)
        axis ([t(1) t(end) yLim(1) yLim(end)]);
        if roi_selec==1
            if jj==1
                c='frontal ROI';
            elseif jj==2
                c='left ROI';
            elseif jj==3
                c='right ROI';
            else
                c='rear ROI';
            end
        elseif roi_selec==2
            c=namroi{jj};
        end
        ylabel(c)
        grid on
    end
    hold on;
    if ~isempty(sss)
    
    hold on;
    px=[t(negmax-150)  t(negmax-150)  t(negmax+200)  t(negmax+200)];
    py=[def.swsd.dispscale -def.swsd.dispscale -def.swsd.dispscale def.swsd.dispscale];
    col=[0.9 0.9 0.1];
    p=patch(px,py,col,'FaceAlpha', 0.5);
    hold on;
    plot(t,0,'k','LineWidth',2)
    end
end
if uiscroll == 1
    uiscroll([1 size(dispdata,2)],10,[],[],yLim);
end


%--------------------------------------------------------------------------
%-----------------    SUBFUNCTION TO MAP THE DELAYS   ---------------------
%--------------------------------------------------------------------------
    
function pm_map(SW,data,origin_count,choice,value,fig,dim, el_set_sphuni,def)


% choice 1 shows delays on a 'hot' colour map of the scalp. the first
% electrode is a pink star. Gray represents electrodes that don't detect the
% considered wave.
% choice 2 draws the data on the scalp in a jet colormap of potentials and
% traces the main trajectory of the wave.

data_channels_eeg=find(strcmpi('eeg',{data.channels(:).type}));
all_names={data.channels(:).label};
channels_eeg_labels=all_names(data_channels_eeg)';
el_label=channels_eeg_labels;
dr=3;

switch choice
    case 1
        
        ed=ones(size(el_label,1),1);
        sv2= zeros(size(el_label,1),1);
        for ielec = 1:size(SW(value).electrodes,2)
            for ichannel = 1:size(el_label,1)
                if strcmpi(char(SW(value).electrodes(ielec)), ...
                            char(el_label{ichannel}))
                   ed(ichannel,1)=SW(value).delays(ielec)-SW(value).delays(1)+dr;
                   sv2(ichannel,1)= 1;
                end
            end
        end
        if dim==0
        %uses interpolation on the unit sphere

       

        SV = crc_CreateSV(ed,[],0,...
                   el_set_sphuni',1,7*pi/12,40,el_label);
        SV2= crc_CreateSV(sv2,[],0,...
                   el_set_sphuni',1,7*pi/12,40,el_label);
        mask = abs(SV2.interpolated)>0.85;
        SV.interpolated=SV.interpolated .* mask;
        label_fig=[char(SW(value).electrodes(1)), ' , type ', data.trials.events(value).type];
        crc_DrawSV(SV,1,0,{label_fig},0,0,6.5); 
        colormap(hot)
        caxis([0 max(ed)]);
        spm_figure('ColorMap','Invert')
        axis off
        cc=colormap;
        cc(1,:)=[0.5 0.5 0.5];
        colormap(cc)
        axis off
        A=get(fig,'Children');
        set(A(2),'Visible','off')
        set(A(3),'Visible','off')
        
        else
        % uses interpolation on a unit disk
        load CRC_electrodes.mat % giopia
        el_names_temp=cell(size(names,2),1);
        for i=1:size(names,2)
            el_names_temp{i}=upper([names{i}]);
        end 
        for i=1:size(origin_count,1)
            origin_count(i,1)=upper([origin_count(i,1)]);
        end 
        [dumb1,dumb2,index2]=intersect(origin_count(:,1),el_names_temp);
        eeg_chan=index2(crc_types(index2)>-2);
        pos_eeg_chan=pos(eeg_chan,:);  
        interp2D(pos_eeg_chan,ed(dumb2),sv2(dumb2),1);
        xlim([0,1])
        ylim([0,1])
        end
          
          
    
    case 2
        if dim==0
        %uses spline interpolation on the unit sphere
          SV = crc_CreateSV(data.data.y(data_channels_eeg,SW(value).start-50:SW(value).start+250),20:220,0,...
                    el_set_sphuni',1,pi/2,10,el_label);
         crc_DrawSV(SV,1,0,SW(value).electrodes(1),1,0,6.5); 
         axis off
        else
            disp('Option not available yet, use delay maps')
        end

end

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INITIALIZE COUNT ON CHANNELS   ---------------
%--------------------------------------------------------------------------
function [origin_count] = pm_origin_count(data)

data_channels_eeg=find(strcmpi('eeg',{data.channels(:).type}));
all_names={data.channels(:).label};
origin_count=all_names(data_channels_eeg)';
for i=1:size(data_channels_eeg,2)
origin_count(i,2)={0};
origin_count(i,3)={0};
end
origin_count(:,2)={0};
origin_count(:,3)={0};

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INTERPOLATE in 2D   --------------------------
%--------------------------------------------------------------------------

function zv=interp2D(el_set_2d,delval,maskval,dra,dst)

% Script for 2D interpolation using matlab built in functions
% Written by Y. Leclercq & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium

if nargin<5
    dst=0.02;
end

% Create space
cc = [0 0];
X=el_set_2d(:,1);
Y=el_set_2d(:,2);
R1 = max(X);
R2=max(Y);
if R1>=R2
    [XI,YI] = meshgrid(-R1:dst:R1);
else
    [XI,YI] = meshgrid(-R2:dst:R2);
end

V=delval;
M=maskval;

% Interpolation & display
[xl,yl,zl] = griddata(X,Y,V,XI,YI,'linear');
[xv,yv,zv] = griddata(X,Y,V,XI,YI,'v4');
[xma,yma,zma]=griddata(X,Y,M,XI,YI,'cubic');
zm = zl; zm(~isnan(zl)) = 1;
zv=zv.*zm;
zmask=zma;
zmask(zma<0.8)=nan;
zv=zv.*zmask;
if dra==1
    pcolor(xv,yv,zv)
    axis image
    hold on
    plot(X,Y,'x')
    colormap(jet)
    spm_figure('Colormap','Invert')
    cc=colormap;
    cc(1,:)=[0 0 0];
    colormap(cc)
    colorbar
    shading interp 
    title('Interpolation using spline, mask within electrodes. Starting area of wave in black')
else
end
return

%_________________________________________________________________________________
function [ xf, yf ] = CartToFlat( x, y, z )
% Convert Cartesian coordinates on surface of unit sphere to 2D Cartesian 
% Flat Map coordinates


rh = sqrt( x.*x + y.*y );

zero = find( rh == 0 );
rh( zero ) = 1e-20;

costheta = x ./ rh;
sintheta = y ./ rh;
    
% Modification by chrisp@fil, 02/02/24,
% to take into account electrodes with phi<-90 or phi>90
%---------------------------------------------------------------------
lz = find(z<0) ;
phi   = atan2( sqrt( x.*x + y.*y ) , z );
% phi(lz) = pi + atan(  sqrt( x(lz).*x(lz) + y(lz).*y(lz) )  ./  z(lz)  );

xf = costheta .* phi;
yf = sintheta .* phi;

