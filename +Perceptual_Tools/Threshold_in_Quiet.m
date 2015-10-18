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
    frequency = frequency' / 1000;
    SPL = 3.64 * frequency.^(-0.8) - 6.5 * exp( -0.6 * (frequency-3.3).^2 ) + 10^-3 * frequency.^4;
    
    
elseif strcmp(method, 'ISO226')
    %% International Organization for Standardization:
    % [2]	B. ISO, “226: 2003:‘Acoustics—Normal equalloudness-level contours,’” Int. Organ. Stand., 2003.
    [spl, F] = iso_226(0);
    
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


function [spl, freq] = iso_226(phon)

f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
     1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];

af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
      0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
      0.243 0.242 0.242 0.245 0.254 0.271 0.301];

Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
       -2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
        2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1];

Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
       11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
       -6.0  -5.4  -1.5   6.0  12.6  13.9  12.3];

%Error Handling
if((phon < 0) || (phon > 90))
    disp('Phon value out of bounds!')
    spl = 0;
    freq = 0;
else
    %Setup user-defined values for equation
    Ln = phon;

    %Deriving sound pressure level from loudness level (iso226 sect 4.1)
    Af=4.47E-3 * (10.^(0.025*Ln) - 1.15) + (0.4*10.^(((Tf+Lu)/10)-9 )).^af;
    Lp=((10./af).*log10(Af)) - Lu + 94;

    %Return user data
    spl = Lp;  
    freq = f;
end
end

