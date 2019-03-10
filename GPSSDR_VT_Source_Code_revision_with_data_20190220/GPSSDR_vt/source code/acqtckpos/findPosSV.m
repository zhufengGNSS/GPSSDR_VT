function  posSV = findPosSV(file,Acquired,ephemeris)
%Purpose
%    
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	Acquired 	- acquisition results
%	ephemeris         - ephemeris 
%Outputs:
%	posSV     	- satellites that can be used to calculation user position
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
sv = Acquired.sv;
posSV           = [];
maxEphUptTime   = 0;
idx             = 0;
for svindex = 1 : length(Acquired.sv)
    prn = sv(svindex);
    if ephemeris(prn).updateflag == 1
        idx = idx + 1;
        posSV = [posSV prn];
        if ephemeris(prn).updatetime(1) > maxEphUptTime
            maxEphUptTime = ephemeris(prn).updatetime(1);
        end
        nAcquired.sv(idx)           = prn;
        nAcquired.SNR(idx)          = Acquired.SNR(svindex);
        nAcquired.Doppler(idx)      = Acquired.Doppler(svindex);
        nAcquired.codedelay(idx)    = Acquired.codedelay(svindex);
        nAcquired.fineFreq(idx)     = Acquired.fineFreq(svindex);
    end
end
save(['nAcquired_',file.fileName],'nAcquired');