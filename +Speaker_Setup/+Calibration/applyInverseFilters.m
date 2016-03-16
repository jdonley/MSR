function [ irs ] = applyInverseFilters( Signals, invIRs )
%APPLYCOMPACTINGFILTER 

irs = Tools.fconv( Signals, invIRs );

end

