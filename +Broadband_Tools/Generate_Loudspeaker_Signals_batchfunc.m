function Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Speech_Setup, Masker_Setup)

%% Initialise
tic;
C = clock;
fprintf('Started execution at %.0f:%.0f:%.0f on the %.0f/%.0f/%.0f\n',C([4:6 3:-1:1]))


%% Analyse broadband signals script
if nargin < 5
    Leakage_Angle = 15;
end

%LUT_res = '512f_256w';
LUT_res = '512f_32w';


Input_file_path = '+Miscellaneous\+Speech_Files\';
%Input_file_path = '+Miscellaneous\+Speech_File_Test\';

%Input_file_path = '+Miscellaneous\+TestAudio_Files\';
%Input_file_path = '+Miscellaneous\+STIPA_Test\';
%Input_file_path = '+Miscellaneous\+Impulse_Response\';
%Input_file_path = '+Miscellaneous\+Sine_Sweep\';

files = Tools.getAllFiles(Input_file_path);


%%
fprintf('\n====== Building Broadband Multizone-Soundfield Loudspeaker-Signals ======\n');
fprintf(['         Zone Weights: ' strrep(sprintf(strrep(repmat('%d',1,length(Zone_Weights)),'d%','d %'),Zone_Weights),' ',' & ') '\n']);
fprintf(['    Noise Mask Levels: ' [strrep(sprintf(strrep(repmat('%d',1,length(Noise_Mask_Levels)),'d%','d %'),Noise_Mask_Levels),' ','dB & ') 'dB'] '\n']);
fprintf(['          Masker Type: ' Masker_Type '\n']);
fprintf(['      Planewave Angle: ' num2str(Speech_Setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle) '\n\n']);n=0;h=[];
fprintf('\tCompletion: ');

F=length(files);
for file = 1:F
    try
        audioread(files{file});
    catch err
        if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
            continue; % Skip unsupported files
        end
    end
    
    W = length(Zone_Weights);
    for w = 1:W
        
        if strcmp(Masker_Type, 'NoMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT( ...
                    files{file}, ...
                    LUT_res, ...
                    Zone_Weights(w), ...
                    Speech_Setup);
        end
        M = length(Noise_Mask_Levels);
        for m = 1:M
            
            if strcmp(Masker_Type, 'FlatMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_FlatMask( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    Speech_Setup);
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask') 
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    Speech_Setup); 
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    Speech_Setup, ...
                    Masker_Setup); 
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl_OffsetNoise')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl_OffsetNoise( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    Speech_Setup, ...
                    Masker_Setup);     %Leakage Angle into Quiet Zone
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl_StereoNoise')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl_StereoNoise( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    Speech_Setup, ...
                    Masker_Setup);     %Leakage Angle into Quiet Zone
                
            end
            
            [n,h] = Tools.showTimeToCompletion( ((file-1)*W*M + (w-1)*M + m) / (F*W*M), n, h);
        end
        
    end
    
    
end


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');