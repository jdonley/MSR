function [N, Frequencies] = Orthogonal_Planewave_Selection( Num_Freqs, N_min, N_max, f_min, f_max )

Frequencies = logspace(  log10(f_min), ...
                         log10(f_max), ...
                         Num_Freqs) ;
                     
N = ceil((N_max - N_min)*(f_max - f_min)^(-1.2)*(Frequencies - f_min).^1.2 + N_min);

end

