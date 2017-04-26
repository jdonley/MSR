function useMeasuredATF(SYS)
%USEMEASUREDATF Swaps the simulated RIRs (Room ATFs) with the measured ATFs
% 
% Syntax:	[OUTPUTARGS] = USEMEASUREDATF(INPUTARGS) Explain usage here
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
% Copyright: Jacob Donley 2017
% Date: 26 April 2017 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    SYS = Current_Systems.loadCurrentSRsystem;
end

%% Find correct setups
ArrType = SYS.system_info.CurrentSpeakerArrayType;
% isCircArray = contains(ArrType, 'circ');
isLineArray = contains(ArrType, 'line');
I = cellfun(@(x) contains(x,ArrType) ,{SYS.Main_Setup.Speaker_Array_Type},'un',0);
S = SYS.Main_Setup([I{:}]);


%% Load current RIRs
[DB,~,RIR_FilePath] = Room_Acoustics.loadRIRDatabaseFromSetup( S(1), SYS.Room_Setup, SYS.system_info.Drive );
if ~(isfield(DB,'isRecordedRIR') && ~DB.isRecordedRIR)
    
    % Get path strings
    [fpath,fname,fext]=fileparts(RIR_FilePath);
    RIRfiles = Tools.getAllFiles(fpath);
    DBs = cellfun(@load , RIRfiles, 'un', 0);
    
    hasIsRecField = cellfun(@(d) isfield(d,'isRecordedRIR'), DBs,'un',0);
    DBsWithIsRecFld = [DBs{[hasIsRecField{:}]}];
    SimRIRs = ~[DBsWithIsRecFld.isRecordedRIR];
    
    if any(SimRIRs)
        DB = DBsWithIsRecFld(SimRIRs);
    end
end
RIRs = DB.RIRs;


%% Get Saved Transfer Functions
filter_location = Tools.getAllFiles([SYS.system_info.Drive SYS.system_info.FilterData_dir]);
filter_location(~contains( filter_location, 'Transfer_Functions'))=[];
filter_location = sort(filter_location);
filts = load( filter_location{1} ); % 1st element should be the newest (most recent) set of filters
TF = filts.TF;
fs = filts.fs;

%%
rir_sim = RIRs.Bright_RIRs;

if isLineArray
    TFAlign = Speaker_Setup.Calibration.linTFalign;
end

rir_rec=[];
for r=1:size(TF,3)
    tftmp = TF(:,:,r);
    if isLineArray
        tftmp = Tools.fconv( tftmp, TFAlign.');
    end
    for s = 1:size(TF,2)
        tmp = decimate(tftmp(:,s), SYS.system_info.fs/SYS.signal_info.Fs);
        rir_rec(r,:,s) = tmp(1:size(rir_sim,2));
    end
end

rir_rec = repmat(rir_rec,size(rir_sim,1)/size(rir_rec,1),1,1);

%% Rename and overwrite RIR files that exist
newRIRs = RIRs;
newRIRs.Quiet_RIRs = rir_rec;
newRIRs.Bright_RIRs = flip(rir_rec,3);

% Re-save simulated RIRs
isRecordedRIR = false;
save([fpath filesep 'Simulated_' fname fext],'RIRs','isRecordedRIR');

% Save new recorded RIRs
RIRs = newRIRs; isRecordedRIR = true;
save(RIR_FilePath,'RIRs','isRecordedRIR');
save([fpath filesep 'Recorded_' fname fext],'RIRs','isRecordedRIR');


end
