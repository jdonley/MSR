function [ match ] = checkMatchingReceivers( setup, room, database_workingdir )
%CHECKMATCHINGRECEIVERS Checks to see if RIR database receiver positions
% match up with saved database.
%
% Returns: TRUE (1) if databases match and FALSE (0) if they do not match.
%

DB1 = Room_Acoustics.loadRIRDatabaseFromSetup( setup, room, database_workingdir );
if isfield(DB1.RIRs,'Matched_Receivers')
    DB2 = load(DB1.RIRs.Matched_Receivers);
else
    error('Database does not contain a path to another database with matching receivers.')
end

match = all(DB1.RIRs.Bright_Receiver_Positions(:) == DB2.RIRs.Bright_Receiver_Positions(:)) ...
     && all(DB1.RIRs.Quiet_Receiver_Positions(:)  == DB2.RIRs.Quiet_Receiver_Positions(:) );

end

