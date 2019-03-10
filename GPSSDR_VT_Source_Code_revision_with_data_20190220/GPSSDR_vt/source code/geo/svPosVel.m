function [sv_xyz, sv_vel, clkcorr_m, clkcorr_m_vel, grpdel] = svPosVel(prn,ephemeris,t,eph_idx)
%Purpose
%   calculate SV position, SV velocity, etc.
%Inputs: 
%	prn         - PRN number
%	ephemeris	- ephemeris
%	t           - time of transmition
%	eph_idx   	- ephemeris index
%Outputs:
%	sv_xyz     	- satellite position in ECEF coordinate
%	sv_vel   	- satellite velocity in ECEF coordinate
%	clkcorr_m  	- satellite clock bias correction in meters
%	clkcorr_m_vel 	- satellite clock drift correction in meters per second
%	grpdel          - group delay in seconds
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
% changing the name of parameters to fit the satnav toolbox
SQRTSMA = ephemeris(prn).sqrta(eph_idx);
DELTAN = ephemeris(prn).deltan(eph_idx);
TOE = ephemeris(prn).toe(eph_idx);
MZERO = ephemeris(prn).M0(eph_idx);
ECCEN = ephemeris(prn).ecc(eph_idx);
ARGPERI = ephemeris(prn).w(eph_idx);
CUS = ephemeris(prn).Cus(eph_idx);
CUC = ephemeris(prn).Cuc(eph_idx);
CRS = ephemeris(prn).Crs(eph_idx);
CRC = ephemeris(prn).Crc(eph_idx);
CIS = ephemeris(prn).Cis(eph_idx);
CIC = ephemeris(prn).Cic(eph_idx);
IZERO = ephemeris(prn).i0(eph_idx);
IDOT = ephemeris(prn).idot(eph_idx);
OMEGAZERO = ephemeris(prn).omegae(eph_idx);
OMEGADOT = ephemeris(prn).omegadot(eph_idx);

toc = ephemeris(prn).toc(eph_idx);
AF0 = ephemeris(prn).af0(eph_idx);
AF1 = ephemeris(prn).af1(eph_idx);
AF2 = ephemeris(prn).af2(eph_idx);
TGD = ephemeris(prn).TGD(eph_idx);


tkc = t - toc;
iter=0;
while tkc > 302400 
    tkc = tkc - 604800;
    iter=iter+1;
    if iter > 3, error('Input time should be time of week in seconds'), end
end
iter=0;
while tkc < -302400 
    tkc = tkc + 604800;
    iter=iter+1;
    if iter > 3, error('Input time should be time of week in seconds'), end
end

    F = -4.442807633e-10;
     
    clkcorr = (AF0 + AF1*tkc + AF2*tkc*tkc) -TGD;  

    
% code of svposeph
gpsPi = 3.1415926535898;
mu = 3986005e8;
OMGedot = 7.2921151467e-5;

    tk = (t-clkcorr) - TOE;
      iter=0;
      while tk > 302400
         tk = tk - 604800;
         iter=iter+1;
         if iter > 3, error('1  Input time should be time of week in seconds'), end
      end
      iter=0;
      while tk < -302400
         tk = tk + 604800;
         iter=iter+1;
         if iter > 3, error('2  Input time should be time of week in seconds'), end
      end
      
    A = ( SQRTSMA )^2;
	n_o = sqrt(mu/(A*A*A));
    n = n_o + DELTAN;      

   Mk = MZERO + n*tk;
   Mk   = rem(Mk + 2*gpsPi, 2*gpsPi);   
   Ek = Mk; sep = 1; oldEk = Ek;
   iter = 0;
   while sep > 1e-13
      Ek = Mk + ECCEN*sin(Ek);
      sep = abs(Ek - oldEk);
      oldEk = Ek;
      iter = iter + 1;
      if iter > 10, break, end
   end
   Ek   = rem(Ek + 2*gpsPi, 2*gpsPi);   
   cos_Ek = cos(Ek);
   sin_Ek = sin(Ek);   
   c1 = 1 - ECCEN * cos_Ek;
   Ek_dot = n/c1;
   
   c2 = sqrt(1-ECCEN*ECCEN);
   sin_vk = ( c2*sin_Ek )/( 1 - ECCEN*cos_Ek );
   cos_vk = ( cos_Ek - ECCEN )/( 1 - ECCEN*cos_Ek );
   vk = atan2(sin_vk,cos_vk);
   vk_dot = Ek_dot * c2/c1;
   
   PHIk = vk + ARGPERI; % remember ARGPERI's unit is rad
   PHIk = rem(PHIk, 2*gpsPi);
   
   c2phik = cos(2*PHIk);
   s2phik = sin(2*PHIk);
   
   delta_uk = CUS*s2phik + CUC*c2phik;
   delta_rk = CRS*s2phik + CRC*c2phik;
   delta_ik = CIS*s2phik + CIC*c2phik;
   
   uk = PHIk + delta_uk;
   uk_dot = vk_dot*(1+2*((CUS*c2phik - CUC*s2phik)));
   
   rk = A*(1-ECCEN*cos_Ek) + delta_rk;
   rk_dot = A*ECCEN*Ek_dot*sin_Ek + 2*vk_dot*(CRS*c2phik - CRC*s2phik);
   
   ik = IZERO + delta_ik + IDOT*tk;
   ik_dot = IDOT + vk_dot*2*(CIS*c2phik - CIC*s2phik);
   
   cos_uk = cos(uk);
   sin_uk = sin(uk);
   
   xxk = rk*cos_uk;    % satellite position In orbital Plane
   yyk = rk*sin_uk;   
   xxk_dot = rk_dot*cos_uk - uk_dot*rk*sin_uk;
   yyk_dot = rk_dot*sin_uk + uk_dot*rk*cos_uk;
   
    OMGk = OMEGAZERO + (OMEGADOT-OMGedot)*(tk) - OMGedot*TOE;
    OMEGADOT = OMEGADOT-OMGedot;
    OMGk = rem(OMGk + 2*gpsPi, 2*gpsPi);
    
	cosO = cos(OMGk);
	sinO = sin(OMGk);
	cosi = cos(ik);
    sini = sin(ik);
     
	sv_xyz(1) = xxk*cosO - yyk*cosi*sinO;
    sv_xyz(2) = xxk*sinO + yyk*cosi*cosO;
    sv_xyz(3) = yyk*sini;
 

    sv_vel(1) = xxk_dot*cosO - OMEGADOT*xxk*sinO - yyk_dot*cosi*sinO + ik_dot*yyk*sini*sinO - OMEGADOT*yyk*cosi*cosO;
    sv_vel(2)= xxk_dot*sinO + OMEGADOT*xxk*cosO + yyk_dot*cosi*cosO - ik_dot*yyk*sini*cosO - OMEGADOT*yyk*cosi*sinO;
    sv_vel(3)= yyk_dot*sini + ik_dot*yyk*cosi;    
    
    c3 = F*ECCEN*SQRTSMA;
    clkcorr_m = 299792458*(AF0 + AF1*tkc + AF2*tkc*tkc + c3*sin_Ek);
    grpdel = TGD;    
    clkcorr_m_vel = 299792458*(AF1 + 2* AF2*tkc + c3*cos_Ek*Ek_dot);