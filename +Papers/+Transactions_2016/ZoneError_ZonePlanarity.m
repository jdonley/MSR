clc;
clear;
close all;
tic;

%%
% SYS = Current_Systems.loadCurrentSRsystem;
SYS = Current_Systems.IEEETransactions_System_H;

%%
f = 1000;
ang = 0;

setup_main = SYS.Main_Setup;
setup_mask = SYS.Masker_Setup;

for s = 1:numel(setup_main)
    
    setup_main(s).Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle
    setup_mask(s).Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle
    
    setup_main(s).Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle = ang;
    setup_mask(s).Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle = ang;
    
    setup_main(s).Multizone_Soundfield.Quiet_Zone = ...
        setup_main(s).Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup_main(s).Multizone_Soundfield.Bright_Zone = ...
        setup_main(s).Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup_main(s).Multizone_Soundfield = setup_main(s).Multizone_Soundfield.setN( -1 ); %Auto set
    setup_main(s).Multizone_Soundfield = setup_main(s).Multizone_Soundfield.createEmptySoundfield('DEBUG');
    if setup_main(s).Loudspeaker_Count>1
        setup_main(s).Multizone_Soundfield = setup_main(s).Multizone_Soundfield.createSoundfield('DEBUG');
    end
    
    
    setup_main(s) = setup_main(s).calc_Loudspeaker_Weights();
    setup_main(s) = setup_main(s).reproduceSoundfield('DEBUG');
end