function [ephemeris] = ini_eph(Acquired)
%Purpose
%   Initialize ephemeris structure
%Inputs:
%	Acquired        - acquisition results
%Outputs:
%	ephemeris       - ephemeris structure 
%---------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
sv = Acquired.sv;

for svindex = 1 : length(sv)
    prn = sv(svindex);
    ephemeris(prn).TOW = [];    
    ephemeris(prn).TOW1 = [];
    ephemeris(prn).sfb1 = [];
    ephemeris(prn).weeknum = [];
    ephemeris(prn).N = [];
    ephemeris(prn).health = [];
    ephemeris(prn).IODC = [];
    ephemeris(prn).TGD = [];
    ephemeris(prn).toc = [];    
    ephemeris(prn).af2 = [];
    ephemeris(prn).af1 = [];
    ephemeris(prn).af0 = [];
    
    ephemeris(prn).IODE2 = [];
    ephemeris(prn).Crs = [];
    ephemeris(prn).deltan = [];
    ephemeris(prn).M0 = [];
    ephemeris(prn).Cuc = [];
    ephemeris(prn).ecc = [];
    ephemeris(prn).Cus = [];
    ephemeris(prn).sqrta = [];
    ephemeris(prn).toe = []; 

    ephemeris(prn).Cic = [];
    ephemeris(prn).omegae = [];
    ephemeris(prn).Cis = [];
    ephemeris(prn).i0 = [];
    ephemeris(prn).Crc = [];
    ephemeris(prn).w = [];
    ephemeris(prn).omegadot = [];
    ephemeris(prn).IODE3 = [];
    ephemeris(prn).idot = [];
    
    ephemeris(prn).updatetime = []; % unit:  ms
    ephemeris(prn).updatetime_tow = []; % unit: second (+6 because of it's next frame)
    
end % end for svindex
end % end function