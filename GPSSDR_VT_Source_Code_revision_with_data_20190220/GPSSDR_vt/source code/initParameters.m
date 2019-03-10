function [file, signal, acq, track, solu, cmn] = initParameters()
%Purpose:
%   Parameter initialization
%Inputs: 
%	None
%Outputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%	acq         - parameters related to signal acquisition,a structure
%	track       - parameters related to signal tracking,a structure
%	solu        - parameters related to navigation solution,a structure
%	cmn         - parameters commmonly used,a structure
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% File parameters
file.fileName       = 'data_20180930_KAITOK_dynamic_f'; 
file.fileRoute      = ['D:\PolyU\GNSS raw data\nslstereo raw data\20180930_KAITOK\',file.fileName,'.dat']; % Change the route of the raw IF data 
file.skip        	= 5000;   	
file.fid           	 = fopen(file.fileRoute,'r','ieee-be');
file.skiptimeVT     = 1000; % when vector tracking begins, unit: ms
file.dataType       = 1;    %1:I; 2:I/Q
file.dataPrecision  = 1;    %1:int8; 2; int16 


%% Signal parameters
signal.IF               = 6.5e6; % unit: Hz 
signal.Fs               = 26e6;	
signal.Fc               = 1575.42e6; 	
signal.codeFreqBasis	= 1.023e6;  	
signal.ms               = 1e-3; % unit: s
signal.Sample           = ceil(signal.Fs*signal.ms);	
signal.codelength       = signal.codeFreqBasis * signal.ms;


%% Acquisition parameters
acq.prnList     = 1:32;     % PRN list
acq.freqStep    = 500;      % unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.L           = 10;       % number of ms to perform FFT


%% Tracking parameters
track.CorrelatorSpacing  	= 0.5;  % unit: chip
track.DLLBW               	= 2;	% unit: Hz
track.DLLDamp           	= 0.707; 
track.DLLGain            	= 0.1;	
track.PLLBW              	= 20;  	
track.PLLDamp             	= 0.707;
track.PLLGain              	= 0.25; 	
track.msEph                 = 40000; % at least 30 seconds 
track.msToProcessCT       	= 5000; % unit: ms
track.msToProcessVT         = 5000; %
track.pdi                   = 1;


%% Navigation parameters
% initial position of the receiver
solu.iniPos	= [22.314/180 * pi, 114.206/180 * pi, 12.849]; % data_20180930_KAITOK_dynamic_f 
solu.cnslxyz = llh2xyz(solu.iniPos); % initial position in the Cartesian ECEF-frame
solu.rate  	= 1000; % unit: Hz


%% commonly used parameters
cmn.vtEnable  	= 1;            % 0: disable vector tracking; 1:enable vector tracking
cmn.cSpeed      = 299792458;    % speed of light, [m/s]
cmn.doy         = 273;       	% Day of year, 273, 282; 362 184;%


%% ionospheirc model (from rinex)
global ALPHA BETA
ALPHA = [0.1024E-07  0.7451E-08 -0.5960E-07 -0.5960E-07]; % 2018/09/30 KAIT0K,HONGKONG
BETA  = [0.8806E+05  0.0000E+00 -0.1966E+06 -0.6554E+05];
 
  