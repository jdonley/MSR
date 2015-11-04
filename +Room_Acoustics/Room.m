classdef Room
    %ROOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Room_Size = [10 10 10];
        Room_Size_txt;
        
        Reproduction_Centre;
        Reproduction_Centre_txt;
        
        Room_Dimensions = 3;
        
        Wall_Absorb_Coeff = 1.0;
        NoReceivers = 32;
    end
    
    methods
        function obj = Room()
            obj = obj.setRoomSize( obj.Room_Size );
            obj = obj.setReproductionCentre( obj.Room_Size ./ 2 );
        end
        
        function obj = setRoomSize(obj, room_size)
           obj.Room_Size = room_size;
           obj.Room_Size_txt = strrep(sprintf(strrep(repmat('%d',1,length(room_size)),'d%','d %'),room_size),' ','x');
        end
        
        function obj = setReproductionCentre(obj, reproduction_centre)
            obj.Reproduction_Centre = reproduction_centre;
            obj.Reproduction_Centre_txt = strrep(sprintf(strrep(repmat('%d',1,length(reproduction_centre)),'d%','d %'),reproduction_centre),' ','x');
        end
    end
    
end

