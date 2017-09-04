function [Path, err, Recordings_Path, path_sub_dirs] = getRecordingsPath( SYS_or_setup, database_res, room, signal_info, database_workingdir, method )
% Summary of this function goes here
% 
% Syntax:	[OUTPUTARGS] = TEMPLATE(INPUTARGS) Explain usage here
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 18 August 2017
% Version: 0.2 (18 August 2017)
% Version: 0.1 (22 February 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latest_method = 'new';
SYS_type = 'Current_Systems.SR_System';

%%
if nargin == 1
    if isa(SYS_or_setup,SYS_type)
        SYS = SYS_or_setup;
        signal_info = SYS.signal_info;
        database_res = SYS.system_info.LUT_resolution;
        setup = SYS.Main_Setup;
        room = SYS.Room_Setup;
        database_workingdir = SYS.system_info.Drive;
        if any(size(signal_info.L_noise_mask)~=1)
            signal_info.L_noise_mask = -inf;
            Tools.simpleWarning(...
                ['Different number of noise levels.\n' ...
                 'Returning path for noise level: ' ...
                 num2str(signal_info.L_noise_mask) ' dB\n'] );
        end
        if any(size(signal_info.methods_list)~=1) && isempty(signal_info.method)
            signal_info.method = signal_info.methods_list{signal_info.methods_list_clean(1)};
            Tools.simpleWarning( ...
                ['Different number of methods.\n' ...
                'Returning path for first clean method: ' ...
                signal_info.method '\n'] );
        elseif all(size(signal_info.methods_list)==1) && isempty(signal_info.method)
            signal_info.method = signal_info.methods_list{1};
        end
    else
        error(['First input argument must be of type: ' SYS_type]);
    end
else    
    setup = SYS_or_setup;
    if nargin < 5
        database_workingdir = 'Z:\';
    end
    if nargin < 4
        signal_info = [];
    end
end

if strcmpi(room.SystemType, 'receive')
    PathType = 'Microphone_Signals';
else
    PathType = 'Recordings';
end

if nargin < 6
    method = latest_method;
end

Recordings_Path = [ ...
    database_workingdir ...
    '+' PathType filesep];

%%
err = false;
try
    
    if strcmpi(method, latest_method)
        
        [~,~,~,~,reproduction_info_dirs, spkr_sig_info_dirs] = Broadband_Tools.getSignalPath( setup, signal_info, database_res );
        
        [~,~,room_info_dir1,room_info_dir2] = Room_Acoustics.getRIRDatabasePath( setup, room );
        room_info_dirs = [room_info_dir1, room_info_dir2, filesep];
        
        path_sub_dirs = [reproduction_info_dirs, ...
            room_info_dirs, ...
            spkr_sig_info_dirs];
        
        Path = [Recordings_Path ...
            path_sub_dirs];
        
        
    elseif strcmpi(method, 'old')
        spkr_sig_dir = ['+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc_LUT_' database_res '\'];
        
        Path = [ ...
            Recordings_Path
            '+Reverb__' num2str(room.NoReceivers) 'Rec_' ...
            room.Room_Size_txt 'Dim_' ...
            room.Reproduction_Centre_txt 'Ctr_' ...
            num2str(room.Wall_Absorb_Coeff) 'Ab\' ...
            spkr_sig_dir];
        
    else
        error('Method to load Results path from setup and room is not supported.')
    end
    
catch ex
    switch ex.identifier
        case 'MATLAB:load:couldNotReadFile'
            warning(['Could not load Results path using the ' method ' method.']);
            err = true;
        otherwise
            rethrow(ex)
    end
end


end

