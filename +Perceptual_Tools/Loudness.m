function [ SPL ] = Loudness( frequency, phons )
%Threshold_in_Quiet Returns the sound pressure level in decibels for the
%                   threshold in quiet for a given frequency using ISO226.


%% International Organization for Standardization:
% [1]	B. ISO, “226: 2003:‘Acoustics—Normal equalloudness-level contours,’” Int. Organ. Stand., 2003.
[spl, F] = iso226( phons );

%% Interpolation
temp = repmat(F, [length(frequency) 1]) - repmat(frequency', [1 length(F)]); % center desireable value
temp(temp >= 0) = NaN;
[~, i_low] = min(abs(temp), [], 2);
i_low = i_low';
i_high = zeros(size(i_low));
i_high(i_low ~= length(temp)) = i_low(i_low ~= length(temp)) + 1;

indices(i_low ~= length(temp)) = (frequency - F(i_low(i_low ~= length(temp)))) ./ ...
    (F(i_high(i_low ~= length(temp))) - F(i_low(i_low ~= length(temp)))) + i_low(i_low ~= length(temp));

indices(i_low == length(temp)) = i_low(i_low == length(temp));

SPL = interp1(spl, indices, 'spline')';

end

