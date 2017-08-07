classdef SR_System
% This class provides complete soundfield reproduction system information
% from compatible loudspeaker setup, signal information, hardware system
% information and room information.

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Date: 20 July 2017
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % Objects
        Room_Setup = [];
        Main_Setup = [];
        Masker_Setup = [];
        signal_info = [];
        system_info = [];
        analysis_info = [];
        publication_info = [];
        
    end
    
    properties (Access = private)
        
    end
    
    %% Contructor
    methods
        function obj = SR_System(~)
            
        end
    end
    
    
end

