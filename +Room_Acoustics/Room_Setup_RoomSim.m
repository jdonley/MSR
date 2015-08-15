clear;
clc;
%%

room.dimension              = [ 5 6 4 ];   % room dimension (x,y,z)
room.humidity               = 0.425;         % relative humidity (0,...,1)
room.temperature            = 20;           % room temperature (celsius)


room.surface.frequency      = [];
room.surface.absorption     = repmat(0.3,6,1);
                             

RT60 = estimateRT60(struct('room',room),'Sabine')





%%


room.dimension              = [ 3 4 3 ];   % room dimension (x,y,z)
room.humidity               = 0.425;         % relative humidity (0,...,1)
room.temperature            = 20;           % room temperature (celsius)


room.surface.frequency      = [];
room.surface.absorption     = repmat(0.3,6,1);
                             

RT60 = estimateRT60(struct('room',room),'Sabine')