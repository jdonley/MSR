function [ SPL ] = Threshold_in_Quiet( frequency, method )
%Threshold_in_Quiet Returns the sound pressure level in decibels for the
%                   threshold in quiet for a given frequency using a
%                   particular method.
if nargin < 2
    method = 'approx';
end


if strcmp(method, 'approx')
    %% Approximation of the threshold in quiet:
    % [1]	E. Terhardt, “Calculating virtual pitch,” Hear. Res., vol. 1, no. 2, pp. 155–182, Mar. 1979.
    frequency = frequency / 1000;
    SPL = 3.64 * frequency.^(-0.8) - 6.5 * exp( -0.6 * (frequency-3.3).^2 ) + 10^-3 * frequency.^4;
    
    
elseif strcmp(method, 'ISO226')
    %% International Organization for Standardization:
    % [2]	B. ISO, “226: 2003:‘Acoustics—Normal equalloudness-level contours,’” Int. Organ. Stand., 2003.
    [spl, F] = iso226(0);
    
    frequencies = length(frequency);
    indices = zeros(frequencies,1);
    for f = 1:frequencies
        temp = F - frequency( f );
        temp(temp >= 0) = NaN;
        [~, i_low] = min(abs(temp));
        if i_low ~= length( F )
            i_high = i_low+1;
            indices(f) = (frequency( f )-F(i_low)) / (F(i_high)-F(i_low)) + i_low;
        else
            indices(f) = i_low;
        end
    end
    SPL = interp1(spl, indices, 'spline');
    
else
    error('The method given is not recognised.');
end

end

