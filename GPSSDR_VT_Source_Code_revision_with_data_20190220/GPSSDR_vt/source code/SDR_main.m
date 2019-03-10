%Purpose:
%   Main function of the GPS software-defined receiver (SDR) based on 
%   vector tracking
%
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%%
clc;
clear;
format long g;
addpath geo             % Geo-related functions, e.g. ionospheric correction function
addpath acqtckpos       % Acquisition, tracking, and postiong calculation functions


%% Parameter initialization 
[file, signal, acq, track, solu, cmn] = initParameters();


%% Acquisition 
if ~exist(['Acquired_',file.fileName,'.mat'])
    Acquired = acquisition(file,signal,acq);
    save(['Acquired_',file.fileName,],'Acquired');    
else
    load(['Acquired_',file.fileName,'.mat']);
end
fprintf('Acquisition Completed. \n');

 
%% Do conventional signal tracking and obtain satellites ephemeris
  if  ~exist(['eph_',file.fileName,'.mat']) || ...
          ~exist(['TckResultCT_forEph_',file.fileName,'.mat']) || ...
            ~exist(['sbf_',file.fileName,'.mat'])
    
    % tracking using conventional DLL and PLL
    TckResultCT_forEph = trackingCT(file,signal,track,Acquired);
    
    % navigaion data decode
    [eph, TckResultCT_forEph, sbf] = naviDecode(Acquired, TckResultCT_forEph);
    
    save(['eph_',file.fileName],'eph'); % ephemeris
    save(['sbf_',file.fileName],'sbf'); % first navi bit point and the beginning of subframe 1
    save(['TckResultCT_forEph_',file.fileName,], 'TckResultCT_forEph');
    
else
    load(['TckResultCT_forEph_',file.fileName,'.mat']);
    load(['eph_',file.fileName,'.mat']);
    load(['sbf_',file.fileName,'.mat']);
end


%% Find satellites that can be used to calculate user position
 posSV  = findPosSV(file,Acquired,eph); 
 load(['nAcquired_',file.fileName,'.mat']);
 Acquired = nAcquired;


%% Do positiong in conventional or vector tracking mode
if cmn.vtEnable == 1
    
    % Added by Jordan Krcmaric, so MATLAB does not throw an error if the
    % convential tracking solution has not been run yet.
    [TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn,Acquired, solu.cnslxyz, eph, sbf);
    
    % initilize VT
    load(['eph_',file.fileName,'.mat']);
    load(['sbf_',file.fileName,'.mat']);    
    load (['tckRstCT_',file.fileName,'.mat']);
    load (['navSolCT_',file.fileName,'.mat']);    
    
    [TckResultVT, navSolutionsVT] = trackingVT_POS(file,signal,track,cmn,Acquired,solu.cnslxyz,eph,sbf,TckResultCT,navSolutionsCT);
else
    [TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn,Acquired, solu.cnslxyz, eph, sbf);
end 

fprintf('Tracking and Positioing Completed.\n\n');
 

