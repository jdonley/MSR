function [N, Frequencies] = Orthogonal_Planewave_Selection( Num_Freqs, N_min, N_max, f_min, f_max, spacing )

if nargin < 6
   spacing = 'log'; 
end

if strcmp(spacing,'log')
    Frequencies = logspace( log10(f_min), ...
                            log10(f_max), ...
                            Num_Freqs) ;
elseif strcmp(spacing,'lin')
    Frequencies = linspace( 0, ...
                            f_max, ...
                            Num_Freqs+1) ;
                        % Truncate to frequencies in the range f_low <-> f_high
                        trunc_index_low  = find(Frequencies < f_min , 1, 'last' ) + 1;
                        trunc_index_high = find(Frequencies >= f_max, 1 ) - 1;
                        if isempty(trunc_index_low)
                            trunc_index_low = 1;
                        end
                        if isempty(trunc_index_high)
                            trunc_index_high = length(Frequencies);
                        end
                        
                        Frequencies = Frequencies( trunc_index_low:trunc_index_high );
end

N = ceil((N_max - N_min)*(f_max - f_min)^(-1.2)*(Frequencies - f_min).^1.2 + N_min);

end

