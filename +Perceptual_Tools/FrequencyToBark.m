function [ Bark ] = FrequencyToBark( frequency )
%FREQUENCYTOBARK Returns the Bark Scale number for a given frequency

Bark = 13 * atan(0.76 * frequency/1000) + 3.5 * atan((frequency / 7500).^2);

end

