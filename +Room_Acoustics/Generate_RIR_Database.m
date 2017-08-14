function Generate_RIR_Database( SYS )
% Generates a database of room impulse responses for a Soundfield Reproduction system
%
% Syntax:	Generate_RIR_Database( SYS )
%
% Inputs:
% 	SYS - Soundfield Reproduction system object
%
% Example:
%  	Room_Acoustics.Generate_RIR_Database( ...
%         Current_Systems.loadCurrentSRsystem );
%
% See also: loadCurrentSRsystem

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015-2017
% Date: 15 August 2015
% Revision: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, SYS = Current_Systems.loadCurrentSRsystem; end

%%
RIRs = [];
Setups = [SYS.Main_Setup; ];

DBsetups = 1:length(Setups);
if isfield(SYS.system_info,'DB_indices')
    DBsetups(~reshape([SYS.system_info.DB_indices{1}],1,[]))=[];
end
tic;
for s = DBsetups
    Setup = Setups(s);
%     delete(gcp('nocreate')); % TODO: Check why this is here. Is it necessary to delete the pool each loop?
    
    Setup.Multizone_Soundfield = ...
        Setup.Multizone_Soundfield.createSoundfield('DEBUG');
    Setup = Setup.calc_Loudspeaker_Weights();
    Setup = Setup.reproduceSoundfield('DEBUG');
    
    %% RIR Generation for a particular setup...
    if s == 1 && ~strcmpi( Setup.Loudspeaker_Type, 'parametric')
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos, rec_b, rec_q ] = ...
            Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
            Setup, ...
            SYS.Room_Setup, ...
            SYS.Room_Setup.Wall_Reflect_Coeff, ...
            SYS.signal_info);
        
        RIRs = struct('Bright_RIRs', RIR_B, ...
            'Bright_Receiver_Positions', Rec_Bright_Pos, ...
            'Quiet_RIRs', RIR_Q, ...
            'Quiet_Receiver_Positions', Rec_Quiet_Pos);
        
    elseif s > 1 && ~strcmpi( Setup.Loudspeaker_Type, 'parametric')
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos ] = ...
            Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_rir_generator( ...
            Setup, ...
            SYS.Room_Setup, ...
            SYS.Room_Setup.Wall_Reflect_Coeff, ...
            SYS.signal_info, ...
            RIRs); % Pass in previous RIRs structure to obtain random 
                   % receiver positions that match previous setup
        
        RIRs = struct('Bright_RIRs', RIR_B, ...
            'Bright_Receiver_Positions', Rec_Bright_Pos, ...
            'Quiet_RIRs', RIR_Q, ...
            'Quiet_Receiver_Positions', Rec_Quiet_Pos, ...
            'Matched_Receivers', {RIR_DB_fullpath});
        
        
    elseif strcmpi( Setup.Loudspeaker_Type, 'parametric') 
        % This assumes the parametric is being used as a masker
        % TODO: Determine if setup is a masker correctly
        [RIR_B, RIR_Q, Rec_Bright_Pos, Rec_Quiet_Pos ] = ...
            Room_Acoustics.RIR_Generation.RIR_from_loudspeaker_setup_PALAnechoic( ...
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
    isRecordedRIR = false;
    save( RIR_DB_fullpath, ...
        'RIRs', ...
        'isRecordedRIR');
    
    if isfield(RIRs,'Matched_Receivers')
        DBPath = RIRs.Matched_Receivers;
        DB = load( DBPath );
        DB.RIRs.Matched_Receivers = RIR_DB_fullpath;
        RIRs = DB.RIRs;
        if isfield(DB, 'isRecordedRIR')
            isRecordedRIR = DB.isRecordedRIR;
        end
        save( DBPath, ...
            'RIRs', ...
            'isRecordedRIR');
    end
    
    %%
    %clear
    %
    %     hold on;
    scat1 = scatter(rec_b(:,1),rec_b(:,2),'.g'); hold on
    scat2 = scatter(rec_q(:,1),rec_q(:,2),'.y');
    if isprop(scat1,'MarkMarkerEdgeAlpha') 
        scat1.MarkMarkerEdgeAlpha = 0.2;
    end
    if isprop(scat2,'MarkMarkerEdgeAlpha')
        scat2.MarkMarkerEdgeAlpha = 0.2;
    end
    
    scatter(Rec_Bright_Pos(:,1),Rec_Bright_Pos(:,2),'ob'); hold on;
    scatter(Rec_Quiet_Pos(:,1),Rec_Quiet_Pos(:,2),'^r'); hold on;
    
    src = [];
    [src(:,1), src(:,2)] = ...
        pol2cart( Setup.Loudspeaker_Locations(:,1), ...
        Setup.Loudspeaker_Locations(:,2));
    src = [src ...
        zeros(size(src,1),size(SYS.Room_Setup.Room_Size,2)-2)] ...
        + repmat(SYS.Room_Setup.Reproduction_Centre([2 1 3]), size(src,1),1);
    scatter(src(:,1),src(:,2),'sk');hold off;
    axis equal;
    
    %axis([0 room.Room_Size(1) 0 room.Room_Size(2)]);
    
end
toc
end


