function [ Final ] = Combine_ThreshMask( Threshold, Mask )
%COMBINE_THRESHMASK Combines a threshold and mask to give a final masking curve

     Mask( Mask <= Threshold ) = 0;
Threshold( Mask ~= 0         ) = 0;

Final = Threshold + Mask;

end

