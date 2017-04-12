function [ setup ] = createSetup( settings )
%CREATESETUP Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
p.CaseSensitive = false;

addOptional(p,'resolution',                     100,               @isnumeric); % 1
addOptional(p,'frequency',                      2000,              @isnumeric); % 2
addOptional(p,'reproduction_radius',            1.0,               @isnumeric); % 3
addOptional(p,'numberof_basisplanewaves',       -1,                @isnumeric); % 4
addOptional(p,'bright_weight',                  1.0,               @isnumeric); % 5
addOptional(p,'quiet_weight',                   2.5,               @isnumeric); % 6
addOptional(p,'unattended_weight',              0.05,              @isnumeric); % 7
addOptional(p,'brightzone_radius',              0.3,               @isnumeric); % 8
addOptional(p,'brightzone_pos_angle',           180,               @isnumeric); % 9
addOptional(p,'brightzone_pos_distance',        0.6,               @isnumeric); % 10
addOptional(p,'brightzone_geometry',            'circle',          @ischar);    % 11
addOptional(p,'brightzone_size',                [1 1],             @isnumeric); % 12
addOptional(p,'brightzone_source_angle',        15,                @isnumeric); % 13
addOptional(p,'brightzone_source_dist',         1.0,               @isnumeric); % 14
addOptional(p,'brightzone_source_type',         'pw',              @ischar);    % 15
addOptional(p,'quietzone_radius',               0.3,               @isnumeric); % 16
addOptional(p,'quietzone_pos_angle',            0,                 @isnumeric); % 17
addOptional(p,'quietzone_pos_distance',         0.6,               @isnumeric); % 18
addOptional(p,'quietzone_geometry',             'circle',          @ischar);    % 19
addOptional(p,'quietzone_size',                 [1 1],             @isnumeric); % 20
addOptional(p,'numberof_loudspeakers',          -1,                @isnumeric); % 21
addOptional(p,'loudspeaker_radius',             1.5,               @isnumeric); % 22
addOptional(p,'maximum_frequency',              8000,              @isnumeric); % 23
addOptional(p,'angleto_firstloudspeaker',       90,                @isnumeric); % 24
addOptional(p,'angleof_loudspeakerarc',         180,               @isnumeric); % 25
addOptional(p,'loudspeaker_model',              'Genelec 8010A',   @ischar);    % 26
addOptional(p,'loudspeaker_object',             [] );                           % 27
addOptional(p,'angleof_loudspeakerarrcentre',   180,               @isnumeric); % 28
addOptional(p,'loudspeaker_spacing',            0.01,              @isnumeric); % 29
addOptional(p,'speaker_array_type',             'circle',          @ischar);    % 30
addOptional(p,'room_size',                      [],                @isnumeric); % 31

parse(p, settings{:});

setup = MultizoneSoundfieldSetup( ...       
    p.Results.resolution, ...                   1
    p.Results.frequency, ...                    2
    p.Results.reproduction_radius, ...          3
    p.Results.numberof_basisplanewaves, ...     4
    p.Results.bright_weight, ...                5
    p.Results.quiet_weight, ...                 6
    p.Results.unattended_weight, ...            7
    p.Results.brightzone_radius, ...            8
    p.Results.brightzone_pos_angle, ...         9
    p.Results.brightzone_pos_distance, ...      10
    p.Results.brightzone_geometry, ...          11
    p.Results.brightzone_size, ...              12
    p.Results.brightzone_source_angle, ...      13
    p.Results.brightzone_source_dist, ...       14
    p.Results.brightzone_source_type, ...       15
    p.Results.quietzone_radius, ...             16
    p.Results.quietzone_pos_angle, ...          17
    p.Results.quietzone_pos_distance, ...       18
    p.Results.quietzone_geometry, ...           19
    p.Results.quietzone_size, ...               20
    p.Results.numberof_loudspeakers, ...        21
    p.Results.loudspeaker_radius, ...           22
    p.Results.maximum_frequency, ...            23
    p.Results.angleto_firstloudspeaker, ...     24
    p.Results.angleof_loudspeakerarc, ...       25
    p.Results.loudspeaker_model, ...            26
    p.Results.loudspeaker_object, ...           27
    p.Results.angleof_loudspeakerarrcentre, ... 28
    p.Results.loudspeaker_spacing, ...          29
    p.Results.speaker_array_type, ...           30
    p.Results.room_size);%                      31

end

function setup = MultizoneSoundfieldSetup( ...
    res, ...     1
    f,...        2
    R,...        3
    N,...        4
    Wb,...       5
    Wq,...       6
    Wu,...       7
    rb,...       8
    angb,...     9
    disb,...     10
    zoGeomb,...  11
    zoSizeb,...  12
    srcangb,...  13
    srcdistb,... 14
    srctypeb,... 15
    rq,...       16
    angq,...     17
    disq,...     18
    zoGeomq,...  19
    zoSizeq,...  20
    L,...        21
    Rl,...       22
    fmax,...     23
    phi,...      24
    phiL,...     25
    spkrmod,...  26
    spkrobj,...  27
    phiLcent,... 28
    spkrspace,...29
    arr_type,... 30
    roomSz)%     31
%%
c = 343;
k = (f/c)*2*pi;
M = ceil(k*R);
kmax = (fmax/c)*2*pi;
Mmax = ceil(kmax*R);

%%
quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rq, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rb, srctypeb, 1.0, srcangb, srcdistb);
quiet.res  = res;
bright.res = quiet.res;
quiet.ZoneGeometry  = zoGeomq;
bright.ZoneGeometry = zoGeomb;
quiet  =  quiet.setZoneSize( zoSizeq );
bright = bright.setZoneSize( zoSizeb );
quiet  =  quiet.setDesiredSoundfield(true, 'suppress_output');
bright = bright.setDesiredSoundfield(true, 'suppress_output');

%%
soundfield = Orthogonal_Basis_Expansion.multizone_soundfield_OBE;
soundfield = soundfield.addSpatialZone(quiet,  disq, angq);
soundfield = soundfield.addSpatialZone(bright, disb, angb);
soundfield.Radius = R;

soundfield.BrightZ_Weight     = Wb;
soundfield.QuietZ_Weight      = Wq;
soundfield.UnattendedZ_Weight = Wu;

if N < 1
    N = 2*M + 1;
end
soundfield = soundfield.setN( floor(N) );

%%
setup = Speaker_Setup.loudspeaker_setup;
setup = setup.setRoomSize(roomSz);
setup = setup.addMultizone_Soundfield(soundfield);
setup = setup.setRadius( Rl );

setup.Speaker_Arc_Angle = phiL;
setup.Angle_FirstSpeaker = phi;

if L < 1
    L = ceil( setup.Speaker_Arc_Angle/360 * 2*Mmax + 1 );
end
setup.Loudspeaker_Count = floor(L);

setup = setup.addLoudspeaker_Object(spkrobj);
setup = setup.setLoudspeakerType(spkrmod);
setup.Speaker_Array_Type = arr_type;
setup.Speaker_Array_Centre = phiLcent;
setup.Speaker_Spacing = spkrspace;
if ~isempty(phiLcent) && ~isempty(spkrspace)
    if ~isnan(phiLcent) && ~isnan(spkrspace)
        setup = setup.setLoudspeakerSpacing();
    end
end

setup = setup.calc_Loudspeaker_Locations;
end