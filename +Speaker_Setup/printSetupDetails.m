function [ details ] = printSetupDetails( setup )
% Creates a string to be used with fprintf to display details of the given setup in human readable form.
%   
% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 25 July 2017
% Revision: 0.1
% 

details = ['\n'...
    'Reproduction' '\n' ...
    '\t' 'Resolution:' '\t\t\t' num2str(setup.res) '\n' ...
    '\t' 'Dimensionality:' '\t\t' num2str(setup.Dimensionality) '\n' ...
    '\t' 'Radius:'   '\t\t\t\t' num2str(setup.Multizone_Soundfield.Radius) '\n' ...
    '\t' 'Unattended Weight:' '\t' num2str(setup.Multizone_Soundfield.UnattendedZ_Weight) '\n' ...
    '\n' ...
    'Bright Zone' '\n' ...
    '\t' 'Weight:'   '\t\t\t' num2str(setup.Multizone_Soundfield.BrightZ_Weight) '\n' ...
    '\t' 'Radius:'   '\t\t\t' num2str(setup.Multizone_Soundfield.Bright_Zone.Radius_q) '\n' ...
    '\t' 'Source Angle:' '\t' num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) '\n' ...
    '\t' 'Position:' '\t' 'X = ' num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.X,10)) '\t' 'Y = ' num2str(round(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y,10)) '\n' ...
    '\t' '\t'      '\t\t' 'Angle = ' num2str(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Angle/pi*180) '\t' 'Distance = ' num2str(setup.Multizone_Soundfield.Bright_Zone.Origin_q.Distance) '\n' ...
    '\n' ...
    'Quiet Zone' '\n' ...
    '\t' 'Weight:'   '\t' num2str(setup.Multizone_Soundfield.QuietZ_Weight) '\n' ...
    '\t' 'Radius:'   '\t' num2str(setup.Multizone_Soundfield.Quiet_Zone.Radius_q) '\n' ...
    '\t' 'Position:' '\t' 'X = ' num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X,10)) '\t' 'Y = ' num2str(round(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y,10)) '\n' ...
    '\t' '\t'      '\t\t' 'Angle = ' num2str(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Angle/pi*180) '\t' 'Distance = ' num2str(setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Distance) '\n' ...
    '\n' ...
    'Loudspeakers' '\n' ...
    '\t' 'Model:'      '\t\t\t' setup.Loudspeaker_Type '\n' ...
    '\t' 'Count:'      '\t\t\t' num2str(setup.Loudspeaker_Count) '\n' ...
    '\t' 'Radius:'     '\t\t\t' num2str(setup.Radius) '\n' ...
    '\t' 'Array Type:'   '\t\t' setup.Speaker_Array_Type '\n' ...
    '\t' 'Array Centre:'   '\t' num2str(setup.Speaker_Array_Centre) '\n' ...
    '\t' 'Spacing:'      '\t\t' num2str(setup.Speaker_Spacing) '\n' ...
    '\t' 'Angle to First:' '\t' num2str(setup.Angle_FirstSpeaker) '\n' ...
    '\t' 'Angle of Arc:'   '\t' num2str(setup.Speaker_Arc_Angle) '\n' ...
    '\n'];
end

