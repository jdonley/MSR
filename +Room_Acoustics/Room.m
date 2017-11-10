classdef Room
% This class provides room setup information for soundfield reproductions

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 14 August 2017
% Version: 0.2 (14 August 2017)
% Version: 0.1 (4 November 2015)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties
        Room_Size = [10 10 10];     % [Depth, Width, Height] <=> [ y, x, z ]
        Room_Size_txt;
        Room_Type = '';             % 'Anechoic' or 'Reverberant'
        
        Reproduction_Centre;        % [Depth, Width, Height] <=> [ y, x, z ]
        Reproduction_Centre_txt;
        
        Room_Dimensions = 3;        % 2-Dimensional or 3-Dimensional
        
        Wall_Absorb_Coeff = 1.0;    % [ Ax1, Ax2, Ay1, Ay2, Az1, Az2 ]
        Wall_Absorb_Coeff_txt = ''; % 'Ax1|Ax2|Ay1|Ay2|Az1|Az2'
        Wall_Reflect_Coeff = 0;     % [ Bx1, Bx2, By1, By2, Bz1, Bz2 ]
        Reflection_Order = -1;      % The order of image sources (reflections) to compute. -1 computes the maximum reflection order for the given length of RIR
        NoReceivers = 32;           % Quantity per zone
        ReceiverPositions = [];     % [Width, Depth, Height] <=> [ x, y, z ]
        ReceiverDirectivityPattern = ... % Microphone directivity pattern
                 'omnidirectional'; % 'omnidirectional', 'subcardioid', 
                                    % 'cardioid', 'hypercardioid', 
                                    % 'bidirectional'
        ReceiverOrientations = ...  % The orientation of the microphones
                             [0 0]; % [azimuth elevation] in radians
                                    
        SystemType = '';            % 'Receive' or 'Transmit'
    end
    
    properties (Access = private)
        RStxt_sep = 'x';
        WACtxt_sep = '-';
    end
    
    methods
        function obj = Room()
            obj = obj.setRoomSize( obj.Room_Size );
            obj = obj.setReproductionCentre( obj.Room_Size ./ 2 );
        end
        
        function obj = setRoomSize(obj, room_size)
           obj.Room_Size = room_size;
           obj.Room_Size_txt = strrep(sprintf(strrep(repmat('%g',1,length(room_size)),'g%','g %'),room_size),' ',obj.RStxt_sep);
        end
        
        function obj = setReproductionCentre(obj, reproduction_centre)
            obj.Reproduction_Centre = reproduction_centre;
            obj.Reproduction_Centre_txt = strrep(sprintf(strrep(repmat('%g',1,length(reproduction_centre)),'g%','g %'),reproduction_centre),' ',obj.RStxt_sep);
        end
        
        function obj = setWall_Absorb_Coeff(obj, wall_absorb_coeff)
           obj = obj.setWall_Absorb_Coeff_txt(wall_absorb_coeff);
            if numel(wall_absorb_coeff) == 1
                wall_absorb_coeff = repmat(wall_absorb_coeff,1,6);
            end
           obj.Wall_Absorb_Coeff = wall_absorb_coeff;
           obj.Wall_Reflect_Coeff = sqrt(1-obj.Wall_Absorb_Coeff);
           obj = obj.setRoom_Type();
        end
        
        function obj = setWall_Absorb_Coeff_txt(obj, wall_absorb_coeff)
            AbsStr = num2str(wall_absorb_coeff);
            AbsStr(regexp(AbsStr,'\s*')) = obj.WACtxt_sep;
            obj.Wall_Absorb_Coeff_txt = strrep(AbsStr,' ','');
        end
        
        function obj = setWall_Reflect_Coeff(obj, wall_reflect_coeff)
           obj = obj.setWall_Absorb_Coeff_txt(1 - wall_reflect_coeff.^2);
            if numel(wall_reflect_coeff) == 1
                wall_reflect_coeff = repmat(wall_reflect_coeff,1,6);
            end
           obj.Wall_Reflect_Coeff = wall_reflect_coeff;
           obj.Wall_Absorb_Coeff = 1 - obj.Wall_Reflect_Coeff.^2;
           obj = obj.setRoom_Type();
        end
        
        function obj = setRoom_Type(obj)
           if obj.Wall_Absorb_Coeff >= 1.0
               obj.Room_Type = 'Anechoic';
           else
               obj.Room_Type = 'Reverberant';
           end
        end
        
        function obj = setReceiverPositions(obj, RecPos)
            obj.ReceiverPositions = RecPos;
        end
        
        function obj = setReceiverDirectivity(obj, Directivity)
            if ischar(Directivity)
               switch Directivity(1:2)
                   case 'om'    % omnidirectional
                       Directivity = 'omnidirectional';
                   case 'su'    % subcardioid
                       Directivity = 'subcardioid';
                   case 'ca'    % cardioid
                       Directivity = 'cardioid';
                   case 'hy'    % hypercardioid
                       Directivity = 'hypercardioid';
                   case 'bi'    % bidirectional
                       Directivity = 'bidirectional';
                   otherwise
                       error('Directivity undefined')
               end
            end
            obj.ReceiverDirectivityPattern = Directivity;
        end
    end
    
end

