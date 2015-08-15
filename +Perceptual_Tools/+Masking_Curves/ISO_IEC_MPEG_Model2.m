function [ SPL ] = ISO_IEC_MPEG_Model2( Bark_Maskee, Bark_Masker, Level_Masker )
%ISO_IEC_MPEG_Model2 ISO/IEC MPEG Psychoacoustic Model 2 speading function
%   The ISO/IEC MPEG Psychoacoustic Model 2 spreading function is derived
%   from the Schroeder spreading function and is given by the following:

dz = Bark_Maskee - Bark_Masker;

SPL = Level_Masker + 15.8111389 + 7.5*(1.05*dz + 0.474) ...
        - 17.5*(1.0 + (1.05*dz + 0.474).^2).^0.5 ...
            + 8*min(0,(1.05*dz - 0.5).^2 - 2*(1.05*dz - 0.5));

end

