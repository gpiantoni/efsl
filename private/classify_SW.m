function [SWtype, SWparam] = classify_SW(cfg, D)
% [SWtype, SWparam] = classify_SW(cfg, D)
%
% where cfg has:
%  .SWest = 'swstream2' (how to estimate streams)
%  .SWcla = 'mostsouth' (how to classify streams)
%
% output is
%  SWtype = a vector (big wave: 1, small wave: 2)
%  SWparam = parametric modulation of the SWtype (f.e. how much south it travels)
%
% extra:
%   you can plot all the streamlines
%   eegpos = coor2D(D, meegchannels(D));
%   figure;
%   plot(eegpos(1,:)', eegpos(2,:)', '.k')
%   hold on
%   classify_SW(cfg, D)

% 11/10/25 removed cfg.compare, added cfg.else
% 11/09/29 mostsouth calculates the relative movement (number of steps), bedinend compares the starting y and ending y (almost always the same, but the second is easier to interpret)
% 11/09/22 include minimum amount of traveling
% 11/09/12 ytype_xylen classifies identical to mostsouth, but parametric is total length
% 11/07/21 moved into savant
% 11/07/18 prevent error if slow waves are not detected
% 11/04/13 added cfg.plotsw
% 11/02/09 no waves: 0, big wave: 1, small wave: 2
% 11/02/08 default with swstream2 and longest_south
% 11/02/07 works
% 11/02/04 created

%---------------------------%
%-defaults

if ~isfield(cfg, 'SWest'); cfg.SWest = 'swstream2'; end
if ~isfield(cfg, 'SWcla'); cfg.SWcla = 'mostsouth'; end
if ~isfield(cfg, 'mintrvl'); cfg.mintrvl = 0; end
if ~isfield(cfg, 'else'); cfg.else = 0; end

% channels that are in common (all apart from ECG and EOG)
if isfield(D.other.CRC.SW, 'DATA4ROI')  % old FASST: only EEG
  warning('using old fasst channel names')

  D_meeg = meeg(D);
  eeglab = chanlabels(D_meeg, meegchannels(D_meeg));
  eegpos = coor2D(D_meeg, meegchannels(D_meeg));
else
  warning('using new fasst channel names')
  not_eeg = find(~ismember({D.channels.type}, {'EEG'}));
  for i = 1:numel(not_eeg)
    D.channels(not_eeg(i)).X_plot2D = NaN;
    D.channels(not_eeg(i)).Y_plot2D = NaN;
  end
  eeglab = {D.channels.label};
  eegpos = [[D.channels.X_plot2D]; [D.channels.Y_plot2D]];
end  
% % simple way of plotting location of slow wave
% [X,Y] = meshgrid(0:.01:1);
% a = griddata( eegpos(1,:), eegpos(2,:), cell2mat(D.other.CRC.SW.origin_count(:,3)), X, Y);
% figure; imagesc(X(1,:), Y(:,1), a)
%---------------------------%

%---------------------------%
%-prepare CONFIG

%-------%
%-cfg for createSWSstream
cfg1 = [];
cfg1.dx = .01;
cfg1.plot = false;
cfg1.plotsw = false;
cfg1.mintrvl = cfg.mintrvl;
cfg1.else = cfg.else;

if isfield(D.other.CRC, 'SW')
  SWtype = zeros(numel(D.other.CRC.SW.SW),1);
else % no slow waves detected
  SWtype = [];
  SWparam = [];
  return
end

progBar = progressBar(numel(D.other.CRC.SW.SW));

for sw = 1:numel(D.other.CRC.SW.SW)
  progBar(sw)
  
  % extra check
  if ~isempty(setxor(eeglab( D.other.CRC.SW.SW(sw).channels), D.other.CRC.SW.SW(sw).electrodes))
    error(['chan labels don''t match chan index in SW ' num2str(sw)])
  end
  
  eegx = eegpos(1, D.other.CRC.SW.SW(sw).channels)';
  eegy = eegpos(2, D.other.CRC.SW.SW(sw).channels)';
  eegz = D.other.CRC.SW.SW(sw).delays;
  
  rmref1 = find(D.other.CRC.SW.SW(sw).channels == 29);
  rmref2 = find(D.other.CRC.SW.SW(sw).channels == 30);
  
  %-----------------%
  % we prob need to remove the two elec that are as ref (not always tho)
  eegx([rmref1 rmref2]) = [];
  eegy([rmref1 rmref2]) = [];
  eegz([rmref1 rmref2]) = [];
  
  %-----------------%
  if cfg1.plot
    crc_SWS_mapping(D,sw,1,[], 1)
  end
  
  if numel(eegx) > 2
    SWstr = feval(cfg.SWest, cfg1, eegx, eegy, eegz);
    
    %-----------------%
    [SWtype(sw), SWparam(sw)] = feval(cfg.SWcla, cfg1, SWstr);
  end
  
end

%-------------------------------------------------------------------------%
%-FUNCTIONS TO ESTIMATE STREAMS
%-------------------------------------------------------------------------%

%---------------------------%
%-swstream2
%---------------------------%
function [SWstr] = swstream2(cfg, eegx, eegy, eegz)
%swstream2
% create some streams using stream2
%
% cfg:
%   cfg.dx = .01;
%
% plot should actually be external

% eegx, eegy: maybe we should remove ref elec that FAST keeps in

%-------%
%-eeg limits: only grid where there are elec
eegxmin = round(min(eegx)/cfg.dx) * cfg.dx;
eegxmax = round(max(eegx)/cfg.dx) * cfg.dx;
eegymin = round(min(eegy)/cfg.dx) * cfg.dx;
eegymax = round(max(eegy)/cfg.dx) * cfg.dx;

[x, y] = meshgrid(eegxmin:cfg.dx:eegxmax, eegymin:cfg.dx:eegymax);

%-------%
%-grid and gradient
V = griddata(eegx, eegy, eegz, x, y);

[u, v] = gradient(V);

%-------%
%-starting points: maybe skip some points
sx = x;
sy = y;

%-----------------%
%-create streamlines
SWstr = stream2(x,y,u,v,sx,sy);

%---------------------------%
%-plot streamline

if cfg.plot
  
  [lenstr] = streamlength(SWstr);
  
  %-------%
  %-which to plot
  toplot = [];
  slenstr = sort(lenstr);
  toplot = lenstr >= slenstr(end-round(numel(slenstr) / 90)); % 90% percentile roughly
  
  %-------%
  %-plot it
  figure
  contour(x,y,V,20)
  hold on
  quiver(x,y,u,v,'k')
  h1 = streamline(SWstr(toplot));
  set(h1,'color','r')
  
end

%-------------------------------------------------------------------------%
%-FUNCTIONS TO CLASSIFY SLOW WAVES
%-------------------------------------------------------------------------%

%---------------------------%
%-mostsouth (if the longest stream goes towards the back)
%---------------------------%
function [SWtype, outlen] = mostsouth(cfg, SWstr)

%-----------------%
%-take longest SW
lenstr = streamlength(SWstr);

[~, maxi] = max(lenstr);
longSW = SWstr{maxi};

%-------%
if cfg.plotsw && strcmp(get(gca, 'nextplot'), 'add') % hold is on
  streamline( {longSW} );
  plot( longSW(1,1), longSW(1,2), '.k', 'markersize', 20)
end

%-------%
%-take away NaN
firstnan = find(isnan(longSW(:,1)), 1);
longSW(firstnan:end,:) = [];

%-------%
%-measure direction
direct = diff(longSW(:,2));
slen = numel(find(direct < 0)); % south-length
nlen = numel(find(direct > 0)); % north-length

if slen - nlen > cfg.mintrvl
  SWtype = 1; % going south
  outlen = slen;
elseif nlen - slen > cfg.mintrvl
  SWtype = 2; % going north
  outlen = nlen;
else
  
  if cfg.else == 1
    SWtype = 1; % going south
    outlen = slen;
  elseif cfg.else == 2
    SWtype = 2; % going north
    outlen = nlen;
  elseif cfg.else == 0
    SWtype = 0;
    outlen = 0;
  end
end

%---------------------------%
%-beginend (compare start and end point: absolute difference)
%---------------------------%
function [SWtype, outlen] = beginend(cfg, SWstr)

%-----------------%
%-take longest SW
lenstr = streamlength(SWstr);

[~, maxi] = max(lenstr);
longSW = SWstr{maxi};

%-------%
if cfg.plotsw && strcmp(get(gca, 'nextplot'), 'add') % hold is on
  streamline( {longSW} );
  plot( longSW(1,1), longSW(1,2), '.k', 'markersize', 20)
end

%-------%
%-take away NaN
firstnan = find(isnan(longSW(:,1)), 1);
longSW(firstnan:end,:) = [];

%-------%
%-measure direction
if longSW(1,2) - longSW(end,2) > cfg.mintrvl
  SWtype = 1; % going south
  outlen = longSW(1,2) - longSW(end,2);
elseif longSW(end,2) - longSW(1,2) > cfg.mintrvl
  SWtype = 2; % going north
  outlen = longSW(end,2) - longSW(1,2);
else
  
  if cfg.else == 1
    SWtype = 1; % going south
    outlen = longSW(1,2) - longSW(end,2);
  elseif cfg.else == 2
    SWtype = 2; % going north
    outlen = longSW(end,2) - longSW(1,2);
  elseif cfg.else == 0
    SWtype = 0;
    outlen = 0;
  end
end

%---------------------------%
%-use y to specify wave type and the whole traveling for parametric
%---------------------------%
function [SWtype, outlen] = ytype_xylen(cfg, SWstr)

%-----------------%
%-take longest SW
lenstr = streamlength(SWstr);

[~, maxi] = max(lenstr);
longSW = SWstr{maxi};

%-------%
if cfg.plotsw && strcmp(get(gca, 'nextplot'), 'add') % hold is on
  streamline( {longSW} );
  plot( longSW(1,1), longSW(1,2), '.k', 'markersize', 20)
end

%-------%
%-take away NaN
firstnan = find(isnan(longSW(:,1)), 1);
longSW(firstnan:end,:) = [];

%-------%
%-measure direction
direct = diff(longSW(:,2));
slen = numel(find(direct < 0)); % south-length
nlen = numel(find(direct > 0)); % north-length

%-------%
%-traveling in all directions
xydiff = diff(longSW); % traveling in small bits
% pythagoras on each small bit
xydist = sqrt(sum(xydiff.^2,2));

xylen = sum(xydist);

if slen - nlen > cfg.mintrvl
  SWtype = 1; % going south
  outlen = xylen;
elseif nlen - slen > cfg.mintrvl
  SWtype = 2; % going north
  outlen = xylen;
else
  
  if cfg.else == 1
    SWtype = 1; % going south
    outlen = xylen;
  elseif cfg.else == 2
    SWtype = 2; % going north
    outlen = xylen;
  elseif cfg.else == 0
    SWtype = 0;
    outlen = 0;
  end
end

%-------------------------------------------------------------------------%
%-EXTRA FUNCTIONS
%-------------------------------------------------------------------------%

%---------------------------%
%-calculate the length of the SW streams
%---------------------------%
function [lenstr] = streamlength(SWstr)

%-------%
%-find longest streams
for st = 1:numel(SWstr)
  firstnan = find(isnan(SWstr{st}(:,1)), 1);
  if isempty(firstnan); firstnan = 0; end
  lenstr(st) = firstnan;
end