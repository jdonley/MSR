function [LoudspeakerSignals_MultiChan, segment_details] = getMSR_Audio( Main_Setup, signal_info, system_info, Masker_Setup, masker_signal_info )
%PLAYMSR Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    Masker_Setup = [];
    masker_signal_info=[];
end


%%
LoudspeakerSignals_MultiChan=[];
segment_details = [];
Setups = [Main_Setup(:); Masker_Setup(:)];
maskers = [zeros(1,numel(Main_Setup)),ones(1,numel(Masker_Setup))];
sig_infos = [signal_info(:); masker_signal_info(:)];
for s = 1:length(Setups)
    if isempty(Setups(s)) || isempty(sig_infos(s))
        break;
    end
    
    %% Get Paths and files
    [Spkr_path,~,~,~,~,~,path_ext] = ...
        Broadband_Tools.getLoudspeakerSignalPath( Setups(s), sig_infos(s), system_info.LUT_resolution, system_info.Drive, 'new');
    
    files = Tools.getAllFiles(Spkr_path);
    files = sort(files);
    
    if ~maskers(s)
        % Continue with only the files that are contained in the original
        % source folder
        files = Tools.keepFilesFromFolder( files, signal_info.speech_filepath);
        if isempty(files), error('No loudspeaker signals found. Have they been generated?'); end
    end
    
    spkr_calib_dir = [system_info.Drive system_info.Calibrated_Signals_dir path_ext];
    
    %% Get signals
    LoudspeakerSignals_MultiChan_ = [];
    fileName_prev = '';
    
    F=length(files);
    for f = 1:F
        
        [~, fileName, fileExt] = fileparts(files{f});
        if ~contains(  fileName, 'Original'  ) && ~signal_info.reference % If not an Original audio file and not a reference recording
           sigType = 'Upsampled';
        elseif contains(  fileName, 'Original'  ) && signal_info.reference
            sigType = 'Original';
        else
            continue;
        end
        
            % Get the file number and file name
            [~,fnameflip] = strtok( flip(fileName), system_info.sc );
            fileName_curr = flip(sscanf(fnameflip,['%*[' system_info.sc ']%s']));
            
            if (isempty(fileName_prev) || ~strcmp( fileName_curr, fileName_prev))
                
               
                y = audioread_extinsensitive([spkr_calib_dir fileName_curr system_info.sc sigType fileExt]);
                if isnan(y)
                    continue;
                end                
                
                % Remove the scaling performed during the calibration
                % which is applied in order to save the audio files without
                % clipping as the signals are changed during the upsampling
                % and filtering process.
                aux_info = load([spkr_calib_dir fileName_curr system_info.sc sigType '.mat']);
                y = y ./ aux_info.scale_value;
                % end rescale
                
                % Remove the scaling performed during the loudspeaker
                % generation process where maskers are scaled to prevent
                % clipping and speech may have also been scaled to prevent
                % clipping if signal_info.inputSignalNorm was set to true.
                common_scaler = 1/db2mag(sig_infos(s).clipLevelAdjust);
                if ~isempty(strfind(sig_infos(s).method, 'Masker'))
                    if sig_infos(s).L_noise_mask <= 0
                        scaler = (db2mag(0)+1); %Plus one is for the amplitude of the clean signal
                    elseif sig_infos(s).L_noise_mask > 0 % For a positive masker we scaled the signals to save up to a 40db noise masker
                        scaler = (db2mag(40)+1);%Plus one is for the amplitude of the clean signal
                    end
                    y = y .* scaler .* common_scaler;
                % The clean input signal may have also been scaled
                elseif sig_infos(s).inputSignalNorm
                    y = y .* common_scaler;
                end
                % end rescale

                
                LoudspeakerSignals_MultiChan_ = [LoudspeakerSignals_MultiChan_; y];
                if s == 1 %speech setup
                    segment_details(end+1).filename = fileName_curr;
                    segment_details(end).fileext = fileExt;
                    segment_details(end).length = size(y,1);
                end
            end
        fileName_prev = fileName_curr;
    end
    if s>1 %This assumes when s>1 we have a masker set of loudspeaker signals (perhaps not a good assumption)
        LoudspeakerSignals_MultiChan_(size(LoudspeakerSignals_MultiChan,1)+1:end,:)=[];
        %LoudspeakerSignals_MultiChan_ = LoudspeakerSignals_MultiChan_ .* db2mag(masker_signal_info.L_noise_mask); % A little bit of a hack %TODO: fix this so speaker signals are leveled correctly        
        % previous line's comment possibly fixed with the following        
%         aux_info = load([spkr_calib_dir fileName_curr system_info.sc sigType '.mat']);
%         LoudspeakerSignals_MultiChan_ = LoudspeakerSignals_MultiChan_ ./ aux_info.scale_value;
    end
    LoudspeakerSignals_MultiChan(:,:,s)=LoudspeakerSignals_MultiChan_;
end
LoudspeakerSignals_MultiChan = sum(LoudspeakerSignals_MultiChan,3);


end


function y = audioread_extinsensitive(file_path)
[p,n,e] = fileparts(file_path);
p=[p,filesep];
if exist([p,n,e],'file')
    audiofile = [p,n,e];
elseif exist([p,n,lower(e)],'file')
    audiofile = [p,n,lower(e)];
elseif exist([p,n,upper(e)],'file')
    audiofile = [p,n,upper(e)];
else
    y = nan;
    return;
end

try
    y = audioread(audiofile);
catch err
    if strcmp(err.identifier, 'MATLAB:audiovideo:audioread:FileTypeNotSupported')
        y=nan; % Skip unsupported files
    end
end
end

