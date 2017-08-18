function [ Path, Name, Ext, err, SubPath, spkr_sig_info_dir, Output_file_path_ext ] = getSignalPath( setup, signal_info, database_res, database_workingdir, SigTypeTxt, method )
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
% Date: 03 February 2016
% Version: 0.1 (03 February 2016)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    method = 'new';
end

if nargin < 4
    database_workingdir = 'Z:\';
    SigTypeTxt = 'IO_Signals';
end

err = false;
try
    %%%%%%%%%%%%%%%%%%%%%%%%% NEW METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
    % Supports database path similarity and better structured path for
    % loudspeaker signal files
    if strcmpi(method, 'new')
        sc = '_'; %Separation Character
        
        Output_file_path     = [database_workingdir '+' SigTypeTxt filesep]; % Can be relative or exact
        
        [~,~,dbsubpath,~,~,~,~,dbfname] = Soundfield_Database.getDatabasePath( setup, database_res );
        [pathstr,fname]=fileparts([dbsubpath,dbfname]);
        SubPath = [pathstr filesep '+' fname filesep];
        spkr_sig_info_dir = [ ...
            '+' num2str(signal_info.f_low ) 'Hz-' ...
            num2str(signal_info.f_high) 'Hz' sc ...
            num2str(signal_info.L_noise_mask) 'dB' sc ...
            num2str(signal_info.weight) 'weight' sc sc ...
            'method' sc signal_info.method filesep ];
        Output_file_path_ext = [SubPath, ...
            spkr_sig_info_dir ];
        Path = [Output_file_path Output_file_path_ext];
        Name     = [signal_info.input_filename sc ];
        Ext      = '.WAV';        
    
    %%%%%%%%%%%%%%%%%%%%%%%%% OLD METHOD %%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmpi(method, 'old')
        
        Output_file_path     = [database_workingdir '+' 'Speaker_Signals' filesep]; % Can be relative or exact
        Output_file_path_ext = ['+' num2str(setup.Radius*2) 'm_SpkrDia\+' num2str(setup.Loudspeaker_Count) 'Spkrs_' num2str(setup.Speaker_Arc_Angle) 'DegArc_LUT_' database_res '\'];        
        Path = [Output_file_path Output_file_path_ext];
        SetupInfo            = ['_' num2str(signal_info.f_low ) 'Hz-' ...
            num2str(signal_info.f_high) 'Hz_' ...
            num2str(setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) 'pwAngle_' ...
            num2str(signal_info.L_noise_mask) 'dB_' ...
            num2str(signal_info.weight) 'weight__with' signal_info.method];
        Name     = [Input_file_name '__' ...
            num2str(setup.Loudspeaker_Count) 'spkrs_' ...
            SetupInfo];
        Ext      = '.WAV';
        
    end
    
catch ex
    switch ex.identifier
        case 'MATLAB:load:couldNotReadFile'
            warning(['Could not get loudspeaker signal path using the ''' method ''' method.']);
            err = true;
        otherwise
            rethrow(ex)
    end
end
end

