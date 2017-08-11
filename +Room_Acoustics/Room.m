classdef Room
    %ROOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Room_Size = [10 10 10];
        Room_Size_txt;
        Room_Type = '';
        
        Reproduction_Centre;
        Reproduction_Centre_txt;
        
        Room_Dimensions = 3;
        
        Wall_Absorb_Coeff = 1.0;
        Wall_Reflect_Coeff = 0;
        NoReceivers = 32; % Per zone
        ReceiverPositions = [];
        
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
        
        function obj = setWall_Absorb_Coeff(obj, wall_absorb_coeff)
           obj.Wall_Absorb_Coeff = wall_absorb_coeff;
           obj.Wall_Reflect_Coeff = sqrt(1-obj.Wall_Absorb_Coeff);
           obj = obj.setRoom_Type();
        end
        
        function obj = setWall_Reflect_Coeff(obj, wall_reflect_coeff)
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
    end
    
end

