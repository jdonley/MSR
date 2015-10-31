function Generate_Loudspeaker_Signals_batchfunc(Zone_Weights, Noise_Mask_Levels, Masker_Type, Planewave_Angle, Leakage_Angle)

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

%loudspeakers = 295;
%loudspeakers = 32;
loudspeakers = 24;
%loudspeakers = 16;

loudspeaker_layout = {'numberof_loudspeakers',        loudspeakers, ...
                      'loudspeaker_radius',           1.5, ...
                      'loudspeaker_model',            'Genelec 8010A', ...
                      'angleof_loudspeakerarrcentre', 180, ...
                      'loudspeaker_spacing',          0.01    };            

zone_layout = {'brightzone_source_angle',     Planewave_Angle};
             
setup = Speaker_Setup.createSetup({...
            zone_layout{:}, ...
            loudspeaker_layout{:}});


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
fprintf(['      Planewave Angle: ' num2str(Planewave_Angle) '\n\n']);n=0;
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
                    setup);
        end
        M = length(Noise_Mask_Levels);
        for m = 1:M
            
            if strcmp(Masker_Type, 'FlatMask')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_FlatMask( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    setup);
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask') 
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    setup); 
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    setup); 
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl_OffsetNoise')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl_OffsetNoise( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    setup, ...
                    Leakage_Angle);     %Leakage Angle into Quiet Zone
                
            elseif strcmp(Masker_Type, 'ZoneWeightMask_AliasCtrl_StereoNoise')
                Broadband_Tools.Loudspeaker_Signal_Calculation.Clean_from_LUT_ZoneWeightedMask_AliasCtrl_StereoNoise( ...
                    files{file}, ...
                    LUT_res, ...
                    Noise_Mask_Levels(m), ...
                    Zone_Weights(w), ...
                    setup, ...
                    Leakage_Angle);     %Leakage Angle into Quiet Zone
                
            end
            
            n = Tools.showTimeToCompletion( ((file-1)*W*M + (w-1)*M + m) / (F*W*M), n);
        end
        
    end
    
    
end


%%
tEnd = toc;
fprintf('\nExecution time: %dmin(s) %fsec(s)\n\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script

% Restart dropbox
%Tools.Dropbox('start');