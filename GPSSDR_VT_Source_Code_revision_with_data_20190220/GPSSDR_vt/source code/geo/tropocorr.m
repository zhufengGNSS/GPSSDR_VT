function tropodel = tropocorr(svxyz,usrxyz)
%TROPOCORR Compute tropospheric correction for
%	the dry component using the Hopfield model
%
%tropodel = tropocorr(svxyz,usrxyz)
%
%	INPUTS
%svxyz = satellite position expressed in ECEF cartesian coordinates
%usrxyz = user position expressed in ECEF cartesian coordinates
%
%	OUTPUTS
%tropodel = estimate of tropospheric error (meters)
 
%Copyright (c) 2002 by GPSoft
%
 
svenu = xyz2enu(svxyz,usrxyz);
el = atan2(svenu(3),norm(svenu(1:2)));
 
tropodel = 2.47/(sin(el)+0.0121);