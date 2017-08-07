function [ TimeSignal ] = OverlapAdd( Frames_x_TimeSig, Overlap )
%OVERLAPADD Summary of this function goes here
%   Detailed explanation goes here

one_minus_overlap_size = floor(size(Frames_x_TimeSig,2)*(1-Overlap));
overlap_size = size(Frames_x_TimeSig,2) - one_minus_overlap_size;

first_few = Frames_x_TimeSig(1, 1:overlap_size);
last_few  = Frames_x_TimeSig(end, one_minus_overlap_size+1:end);

overlapping_1 = Frames_x_TimeSig(1:end-1, one_minus_overlap_size+1:end);
overlapping_2 = Frames_x_TimeSig(2:end  , 1:overlap_size  );

added_part = (overlapping_1 + overlapping_2)';

non_overlapping = Frames_x_TimeSig(:, overlap_size+1:one_minus_overlap_size);

if ~isempty(non_overlapping)
    TimeSignal = [non_overlapping(1:end-1,:) added_part']';
    TimeSignal = TimeSignal(:);
    TimeSignal = [first_few'; TimeSignal; non_overlapping(end,:)'; last_few'];
else
    TimeSignal = [first_few'; added_part(:); last_few'];
end

end

