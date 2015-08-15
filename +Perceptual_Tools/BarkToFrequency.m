function [ frequency ] = BarkToFrequency( Bark )
%FREQUENCYTOBARK Returns the frequency for a given Bark Scale number

frequency = ((( exp( 0.219*Bark ) / 352) + 0.1) .* Bark ...
                - 0.032 * exp( -0.15 * (Bark-5).^2 )) * 1000;

end

