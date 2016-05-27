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
Setups = {Main_Setup; Masker_Setup};
sig_infos = {signal_info; masker_signal_info};
for s = 1:length(Setups)
    if isempty(Setups{s}) || isempty(sig_infos{s})
        break;
    end
    
    %% Get Paths and files
    [Spkr_path,~,~,~,~,~,path_ext] = ...
        Broadband_Tools.getLoudspeakerSignalPath( Setups{s}, sig_infos{s}, system_info.LUT_resolution, system_info.Drive, 'new');
    
    files = Tools.getAllFiles(Spkr_path);
    files = sort(files);
    
    
    spkr_calib_dir = [system_info.Drive system_info.Calibrated_Signals_dir path_ext];
    
    %% Get signals
    LoudspeakerSignals_MultiChan_ = [];
    fileName_prev = '';
    
    F=length(files);
    for f = 1:F
        
        [~, fileName, fileExt] = fileparts(files{f});
        if isempty( strfind( fileName, 'Original' ) ) && ~signal_info.reference % If not an Original audio file and not a reference recording
           sigType = 'Upsampled';
        elseif ~isempty( strfind( fileName, 'Original' ) ) && signal_info.reference
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
                
                LoudspeakerSignals_MultiChan_ = [LoudspeakerSignals_MultiChan_; y];
                if s == 1 %speech setup
                    segment_details(end+1).filename = fileName_curr;
                    segment_details(end).fileext = fileExt;
                    segment_details(end).length = size(y,1);
                end
            end
        fileName_prev = fileName_curr;
    end
    if s>1
        LoudspeakerSignals_MultiChan_(size(LoudspeakerSignals_MultiChan,1)+1:end,:)=[];
        %LoudspeakerSignals_MultiChan_ = LoudspeakerSignals_MultiChan_ .* db2mag(masker_signal_info.L_noise_mask); % A little bit of a hack %TODO: fix this so speaker signals are leveled correctly        
        % previous line's comment possibly fixed with the following        
        aux_info = load([spkr_calib_dir fileName_curr system_info.sc sigType '.mat']);
        LoudspeakerSignals_MultiChan_ = LoudspeakerSignals_MultiChan_ ./ aux_info.scale_value;
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

