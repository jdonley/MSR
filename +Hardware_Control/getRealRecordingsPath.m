function [Path, err] = getRealRecordingsPath( SYS_or_setup, database_res, room, signal_info, database_workingdir, method )
%GETREALRECORDINGSPATH Summary of this function goes here

latest_method = 'new';
SYS_type = 'Current_Systems.SR_System';

%%
if nargin == 1
    if isa(SYS_or_setup,SYS_type)
        wrnCol = [255,100,0]/255;
        SYS = SYS_or_setup;
        signal_info = SYS.signal_info;
        database_res = SYS.system_info.LUT_resolution;
        setup = SYS.Main_Setup;
        room = SYS.Room_Setup;
        database_workingdir = SYS.system_info.Drive;
        if any(size(signal_info.L_noise_mask)~=1)
            signal_info.L_noise_mask = -inf;
            cprintf(wrnCol, ...
                ['Different number of noise levels.\n' ...
                 'Returning path for noise level: ' ...
                 num2str(signal_info.L_noise_mask) ' dB\n'] );
        end
        if any(size(signal_info.methods_list)~=1) && isempty(signal_info.method)
            signal_info.method = signal_info.methods_list{signal_info.methods_list_clean(1)};
            cprintf(wrnCol, ...
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

if nargin < 6
    method = latest_method;
end


RealWorldRecordings_Path = '+Physical_World\';

err = false;
try
    
    if strcmpi(method, latest_method)
        
        [~,~,rec_pre_path,recording_subpath] = Results.getRecordingsPath( ...
            setup, ...
            database_res, ...
            room, ...
            signal_info, ...
            database_workingdir, ...
            method);
        
        Path = [rec_pre_path, RealWorldRecordings_Path, recording_subpath];
                   
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

