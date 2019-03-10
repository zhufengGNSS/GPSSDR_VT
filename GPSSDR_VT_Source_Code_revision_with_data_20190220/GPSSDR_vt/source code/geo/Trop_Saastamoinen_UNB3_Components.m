function [K_dry, K_wet] = Trop_Saastamoinen_UNB3_Components(doy,lat,alt)

UNB3_GM = 9.80665;
UNB3_RD = 287.054;
UNB3_K1 = 0.000077604;
UNB3_K2 = 0.382;

[T, T0, P, WVP, Beta, Lambda] = Get_UNB3_Model(doy, lat, alt);

 K_dry = P * UNB3_K1*UNB3_RD/UNB3_GM;
 K_wet = WVP * UNB3_K2*UNB3_RD /...
                                ((UNB3_GM*(Lambda+1) - Beta*UNB3_RD) * T0);
end