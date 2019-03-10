function [T, T0, P, WVP, beta, lambda] = Get_UNB3_Model(doy,lat,alt)

UNB3_GM = 9.80665;
UNB3_RD = 287.054;
UNB3_K1 = 0.000077604;
UNB3_K2 = 0.382;

R2D = 180/pi;

avg = [15.0 1013.25 299.65 26.31 0.00630 2.77;
       30.0 1017.25 294.15 21.79 0.00605 3.15;
       45.0 1015.75 283.15 11.66 0.00558 2.57;
       60.0 1011.75 272.15  6.78 0.00539 1.81;
       75.0 1013.00 263.65  4.11 0.00453 1.55];

   
amp = [15.0  0.00  0.00 0.00 0.00    0.00;
       30.0 -3.75  7.00 8.85 0.00025 0.33;
       45.0 -2.25 11.00 7.24 0.00032 0.46;
       60.0 -1.75 15.00 5.36 0.00081 0.74;
       75.0 -0.50 14.50 3.39 0.00062 0.30];
   
doy2rad = 2 * pi / 365.25;
ep = UNB3_GM / UNB3_RD;

% lat = lat * R2D; % if degree

% /* Southern hemisphere and yearly variation */
 if(lat < 0.0)
   doy = doy - 211.0;
 else
   doy = doy - 28.0;
 end
 
%  /* phase in annual cycle */
 cosphs = cos(doy*doy2rad);
 
%   /* Initialize pointers to look-up table */
 lat = abs(lat); 
  if(lat >= 75.0)
   p1 = 4;
   p2 = 4;
   m = 0;
  elseif(lat <= 15.0)
   p1 = 0;
   p2 = 0;
   m = 0;
  else
      
   p1 = floor( (lat - 15) / 15) + 1; % +1 by Q-mo
   p2 = p1 + 1;
   
   m = (lat -avg(p1,1)) / (avg(p2,1) - avg(p1,1));
  end
 
  
%    /* avg. surface values */
 Pavg =   m * (avg(p2,2) - avg(p1,2)) + avg(p1,2); 
 Tavg =   m * (avg(p2,3) - avg(p1,3)) + avg(p1,3); 
 WVPavg = m * (avg(p2,4) - avg(p1,4)) + avg(p1,4);
 Bavg =   m * (avg(p2,5) - avg(p1,5)) + avg(p1,5);
 Lavg =   m * (avg(p2,6) - avg(p1,6)) + avg(p1,6);
   
%  /* annual cycle amplitude of surface values*/
 Pamp =   m * (amp(p2,2) - amp(p1,2)) + amp(p1,2); 
 Tamp =   m * (amp(p2,3) - amp(p1,3)) + amp(p1,3); 
 WVPamp = m * (amp(p2,4) - amp(p1,4)) + amp(p1,4);
 Bamp =   m * (amp(p2,5) - amp(p1,5)) + amp(p1,5);
 Lamp =   m * (amp(p2,6) - amp(p1,6)) + amp(p1,6);
   
%  /* total surface values */
 T0   = Tavg - Tamp * cosphs; 
 P0   = Pavg - Pamp * cosphs;
 WVP0 = WVPavg - WVPamp * cosphs;
 beta = Bavg - Bamp * cosphs;
 lambda = Lavg - Lamp * cosphs;

%  /* met values at altitude */
 T = T0 - beta * alt;
 P = P0 * ((T/T0)) ^ ( ep/(beta));
 WVP = WVP0 * ((T/T0)) ^ ((ep * ((lambda) + 1)/(beta))-1);
  
end