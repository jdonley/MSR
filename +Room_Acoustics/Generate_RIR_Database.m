clear;
clc;
% clear classes;
close all;
tic;

%%
SYS = Current_Systems.loadCurrentSRsystem;


RIRs = [];
Setups = [SYS.Main_Setup; ];
for s = 1:length(Setups)
    Setup = Setups(s);
    delete(gcp('nocreate'));
    
    Setup.Multizone_Soundfield = Setup.Multizone_Soundfield.createSoundfield('DEBUG');
    Setup = Setup.calc_Loudspeaker_Weights();
    Setup = Setup.reproduceSoundfield('DEBUG');
    
    %% RIR Generation for a particular setup...
    if ~strcmpi( Setup.Loudspeaker_Type, 'parametric')
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
            Setup, ...
            SYS.Room_Setup, ...
            repmat(SYS.Room_Setup.Wall_Reflect_Coeff,1,6), ...
            SYS.signal_info);
        
        RIRs = struct('Bright_RIRs', RIR_B, ...
            'Bright_Receiver_Positions', Rec_Bright_Pos, ...
            'Quiet_RIRs', RIR_Q, ...
            'Quiet_Receiver_Positions', Rec_Quiet_Pos);
        
    else
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos ] = Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_PALAnechoic( ...
            Setup, ...
            SYS.Room_Setup, ...
            SYS.signal_info, ...
            RIRs, ...
            RIRs.Quiet_Receiver_Positions); %Normalise at the target quiet zone
        
        RIRs = struct('Bright_RIRs', RIR_B, ...
            'Bright_Receiver_Positions', Rec_Bright_Pos, ...
            'Quiet_RIRs', RIR_Q, ...
            'Quiet_Receiver_Positions', Rec_Quiet_Pos, ...
            'Matched_Receivers', {RIR_DB_fullpath});
    end
    
        
    
    %% Save the RIRs to a database for reuse
    RIR_DB_fullpath = Room_Acoustics.getRIRDatabasePath( ...
        Setup, ...
        SYS.Room_Setup, ...
        SYS.system_info.Drive);
    RIRDBpath = fileparts(RIR_DB_fullpath);
    
    if ~exist( RIRDBpath,'dir'); mkdir( RIRDBpath ); end
    save( RIR_DB_fullpath, ...
        'RIRs');
    
    if strcmpi( Setup.Loudspeaker_Type, 'parametric')
        DBPath = RIRs.Matched_Receivers;
        DB = load( DBPath );
        DB.RIRs.Matched_Receivers = RIR_DB_fullpath;
        RIRs = DB.RIRs;
        save( DBPath, 'RIRs');
    end
    
    %%
    %clear
    %
    hold on;
    scatter(rec_b(:,1),rec_b(:,2),'.g'); hold on
    scatter(rec_q(:,1),rec_q(:,2),'.y');
    scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob'); hold on;
    scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'or'); hold on;
    
    src = [];
    [src(:,1), src(:,2)] = pol2cart( Setup.Loudspeaker_Locations(:,1), Setup.Loudspeaker_Locations(:,2));
    src = [src zeros(size(src,1),size(SYS.Room_Setup.Room_Size,2)-2)] + repmat(SYS.Room_Setup.Reproduction_Centre, size(src,1),1);
    scatter(src(:,1),src(:,2),'ok');hold off;
    axis equal;
    %axis([0 room.Room_Size(1) 0 room.Room_Size(2)]);
    
end

%%
toc

    



