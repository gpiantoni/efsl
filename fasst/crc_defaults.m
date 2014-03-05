function ouput=crc_defaults(label)
% Sets the defaults which are used by FAST, fMRI Artefact removal and Sleep
% scoring Toolbox.
%
% FORMAT crc_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure your 'crc_defaults' is the first one found in the 
% path. See matlab documentation for details on setting path.
%
% Care must be taken when modifying this file!
%
% The structure and content of this file are largely inspired by SPM.
%_______________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: crc_defaults.m 210 2010-05-28 08:02:25Z jessica $

switch label

    %%%%%%%%%%%%%%%%%%% Artefact Rejection Parameters %%%%%%%%%%%%%%%%%%%%%%%%%

    % Parameters for the gradient artefact rejection, crc_gar.m
    %-----------------------------------------------

    case 'crc_gar'
        
        crc_defaults.gar.prefix         = 'cga_';
        crc_defaults.gar.output_fs      = 500 ;     % Frequency sampling of the output file
        crc_defaults.gar.Nsc_aver       = 30 ;      % Average computed on ... scans
        crc_defaults.gar.Nsc_skipped_be = 2 ;       % Number of scans skipped before correction
        crc_defaults.gar.UseScanMrk     = 1 ;       % Using volume marker from the scanner (1), or not (0)
        
        % If UseScanMrk = 1;
        crc_defaults.gar.ScanMrk1       = 128 ;     % Marker number coming from the scanner
        crc_defaults.gar.ScanMrk2       = 1002 ;       % Secondary scanner marker.
        crc_defaults.gar.MrkStep        = 1 ;       % Number of marker between 2 successive volumes
                                                    % If you have a scanner marker every slice,
                                                    % use the number of slices !
        
        % If UseScanMrk = 0;
        crc_defaults.gar.TR             = 2.46; % in sec
        crc_defaults.gar.AutoChk        = 1;    % Autodetection of beginning & end
            % If crc_defaults.gar.AutoChk = 1
        crc_defaults.gar.Threshold      = 350;      
        crc_defaults.gar.DetStep        = 1 ;   % Step use for detection (in sec)
        crc_defaults.gar.DetChan        = 1 ;   % Channel use for detection
            % If crc_defaults.gar.AutoChk = 0        
        crc_defaults.gar.beg            = 1 ;       % in sec
        crc_defaults.gar.nd             = 2500;     % in sec

        % Parameters for cICA pulse artefact calculation,
        %-----------------------------------------------
        % used in crc_bcgremoval_rt.m

    case 'crc_par'

        % Filename convention.
        crc_defaults.par.cicapref  = 'cICA_'; %filename prefix added after using cICA
        crc_defaults.par.pcapref   = 'cpca_'; %filename prefix added after using PCA
        crc_defaults.par.aaspref   = 'caas_'; %filename prefix added after using AAS
        crc_defaults.par.iicapref  = 'ciica_'; %filename prefix added after using iICA
        crc_defaults.par.pcaaspref = 'cpcaas_'; %filename prefix added after using PCA/AAS

        % Parameters for pulse artefact rejection
        %----------------------------------------
        crc_defaults.bcgrem.size        = 90 ; % in sec (60 is the minimum)
        crc_defaults.bcgrem.step        = 60 ; % about 2/3 size
        % Channels used to build the constrain vector.
        % They're taken from around the head
        crc_defaults.par.Refnames = ...
            {'AF3' 'AF4' 'FPZ' 'FP1' 'FP2' 'O1' 'O2' 'F1' 'F2' 'C1' 'C2' 'PZ' 'OZ' 'F5' ...
            'F6' 'C5' 'C6' 'P5' 'P6'};
        crc_defaults.par.NitKmeans = 150;

        % Parameters for pulse artefact matrix quality assessment
        %--------------------------------------------------------
        crc_defaults.par.useinitseg=1; % Use initial segment to assess the quality of the correction matrix.
        crc_defaults.par.additioseg=10; % Number of additional random segments.
        %Note that if both of these value equals 0, the software use the
        %1st segment anyway.
        crc_defaults.par.length = 60; % Length of each segment in seconds.
        
        %%%%%%%%%%%%%%%%%%% Display parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Parameters for display one file,
        %-----------------------------------------------
        % used in dis_main.m, dis_cmp.m
    case 'dis_one'
        crc_defaults.one.winsize      = 20; % in sec
        crc_defaults.one.scale        = 150; % in µV
        crc_defaults.one.filtEEG    = [.1 20]; % in Hz
        crc_defaults.one.filtEMG    = [10 125]; % max is adapted to true sampling rate
        crc_defaults.one.filtEOG    = [.1 4]; % in Hz

        % Parameters for display multiple files
        %-----------------------------------------------
    case 'dis_mul'
        
        crc_defaults.mul.winsize=20;
        crc_defaults.mul.scale=150;
        crc_defaults.mul.filtEEG    = [.1 20]; % in Hz
        crc_defaults.mul.filtEMG    = [10 125]; % max is adapted to true sampling rate
        crc_defaults.mul.filtEOG    = [.1 4]; % in Hz
        
        % Parameters for display power spectrum
        %-----------------------------------------------
    case 'dis_frq'
        crc_defaults.disfrq.scale       = [200 0.3 5]; % Absolute Value / Relative Value / Mongrain Value
        crc_defaults.disfrq.subsmpl     = 16;
        crc_defaults.disfrq.maxpix      = 1600;

        % Parameters for scoring file, dis_score.m
        %-----------------------------------------------
    case 'dis_score'
        
        crc_defaults.score.scale        = 150; % in µV

        %%%%%%%%%%%%%%%%%%%% Computing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Parameters for computing spectrogram or power spectrum, crc_spectcompute.m
        %-----------------------------------------------
    case 'cps'
        crc_defaults.cps.uplimit        = 25 ;
        crc_defaults.cps.downlimit      = .5 ;
        crc_defaults.cps.duration       = 4 ;
        crc_defaults.cps.step           = 2 ;
        crc_defaults.cps.scorer         = 1 ;
        crc_defaults.cps.reference      = 1 ;
        crc_defaults.cps.maxmemload     = 50 * 1024^2; % i.e. 50Mb
    
    case 'ds'
        % Parameters for downsampling
        %-----------------------------------------------
        crc_defaults.ds.stepse = 60; % sec
        crc_defaults.ds.fs = 250; % Hz
        crc_defaults.ds.prefix = 'ds_';
        
        % Parameters for Slow Wave detection
        %-----------------------------------------------
    
    case 'swsd'
        crc_defaults.swsd.dispscale=200; %microV
        crc_defaults.swsd.bpfc=200; %lowpass cutoff of bandpass filter
        crc_defaults.swsd.ROI_centers=[0.5 0.6341;  %Fz  %position of ROI centers on the 2D unit disk
                                    0.3536 0.5013; %C3
                                    0.6464 0.5013; %C4
                                    0.5 0.3575]; %Pz
        crc_defaults.swsd.param1=struct('SWlength',[250 1250 1500 200],... % minimum and maximum time durations between down zero crossing and up zero crossing
                'SWmAmpl',[-40 -80 75 140]);                             % Massimini criteria of magnitude for SWS (-80 140) and delta waves (-40 75)
        crc_defautls.swsd.highfc=0.2; %highpass frequency cutoff
        crc_defautls.swsd.lowfc=4;    %lowpass frequency cutoff
        crc_defautls.swsd.stagesw=[3 4]; %stages to extract from original scored file
        crc_defaults.swsd.butterorder=4;  %order of butterworth filter
    otherwise
        warning('Improper use of this function');
end
ouput=crc_defaults;






