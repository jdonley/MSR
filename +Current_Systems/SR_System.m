classdef SR_System
    % SR_SYSTEM - This class provides complete soundfield reproduction
    %             system information from compatible loudspeaker setup,
    %             signal information, hardware system information and room
    %             information.
    %
    %
    %   Author: Jacob Donley, University of Wollongong, Australia
    %   Email: Jacob.Donley089@uowmail.edu.au
    %
    
    properties
        % Objects
        Room_Setup = [];
        Main_Setup = [];
        Masker_Setup = [];
        signal_info = [];
        system_info = [];
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

