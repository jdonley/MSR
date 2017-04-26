function TFalign = linTFalign(SYS)
%LINTFALIGN Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = LINTFALIGN(INPUTARGS) Explain usage here
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 26 April 2017 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    SYS = Current_Systems.loadCurrentSRsystem;
end

Fs = SYS.system_info.fs;

isLineArray = contains(SYS.system_info.CurrentSpeakerArrayType, 'line');
if isLineArray
    I = cellfun(@(x) contains(x,'line') ,{SYS.Main_Setup.Speaker_Array_Type},'un',0);
    if sum([I{:}])~=1, error('Too many or too few line arrays in the given System'); end
    S = SYS.Main_Setup([I{:}]);
    
    r = S.Loudspeaker_Locations(:,2);
    rAlign = S.Radius;
    impi = round((r-rAlign)/S.c*Fs); % Samples to delay
    
    TFalign = zeros( numel(impi), max(impi)-min(impi)+1 );
    
    TFalign( sub2ind(size(TFalign),1:numel(impi),(impi-min(impi)+1).') ) ...
        = abs(exp(1i.*r)./r) ./ abs(exp(1i*rAlign)/rAlign); % Uses 3D acoustic transfer function to adjust power levels
end

end
