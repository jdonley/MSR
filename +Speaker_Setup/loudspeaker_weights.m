function [ speaker_weights ] = loudspeaker_weights( multizone_soundfield )
%LOUDSPEAKER_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here

          phi_0 = 2*pi; 
             Rl = multizone_soundfield.Radius;
             
              M = multizone_soundfield.getGlobalModeLimit;
             Q0 = ceil( phi_0*(2*M+1)/(2*pi) );
speaker_weights = zeros(1, Q0);
              C = zeros(2*M+1,1);
for m = -M:M
       C(m+M+1) = 2 / (1i*pi*besselh(m, multizone_soundfield.k_global * Rl));
end
           Beta =  C .* multizone_soundfield.Alpha_Coeffs;
    delta_phi_s = (phi_0) / Q0; 
          phi_q = phi_0 + ((1:Q0) - 1) * delta_phi_s;
              m = -M:M;

 for q = 1:Q0
  speaker_weights(q) = sum(Beta .* exp(1i*m*phi_q(q)) * delta_phi_s);
 end
end

