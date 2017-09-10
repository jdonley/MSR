function [Time,RIR] = import_RoomImpulseResponse_vs_Time( SYS )
%import_RoomImpulseResponse_vs_Time Import numeric data from a text file as column vectors.
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_SoundPressureLevels_vs_Freq(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up]
%   = import_SoundPressureLevels_vs_Freq(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [BlockLen,SUPP,ConfInt_Low,ConfInt_Up] = import_RoomImpulseResponse_vs_Time('RIR_Results.mat',1, 180);
%


%%

S = load([Results.getResultsPath(SYS) 'RIR_Results.mat']);
RIR = [S.S{1}{1:2}];
Time = (0:size(RIR,1)-1).'/SYS.signal_info.Fs;


        
end


