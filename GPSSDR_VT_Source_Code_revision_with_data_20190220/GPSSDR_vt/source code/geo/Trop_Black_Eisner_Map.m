function [m_dry] = Trop_Black_Eisner_Map(cos_elev)

 m_dry = 1.0/sqrt(1.0 - cos_elev*cos_elev/1.002001);
 
end