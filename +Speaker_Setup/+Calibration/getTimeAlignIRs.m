function [ TimeAlignIRs ] = getTimeAlignIRs( IRs )
%GETTIMEALIGNIRS Summary of this function goes here
%   Detailed explanation goes here

% Time-alignment 
[~,IRdelays] = max(abs(IRs));
IRdelays = max(IRdelays)-IRdelays + 1;

[g1,g2]=meshgrid( IRdelays, 1:max(IRdelays) );
TimeAlignIRs = double((g1-g2)==0);
%TimeAlignIRs = flip(TimeAlignIRs,2);
end

