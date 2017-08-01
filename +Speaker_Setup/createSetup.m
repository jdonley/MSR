function [ setup ] = createSetup( settings )
%CREATESETUP Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
p.CaseSensitive = false;

addOptional(p,'resolution',                     100,               @isnumeric); % 1
addOptional(p,'frequency',                      2000,              @isnumeric); % 2
addOptional(p,'reproduction_radius',            1.0,               @isnumeric); % 3
addOptional(p,'reproduction_geometry',          'circle',          @ischar);    % 4
addOptional(p,'reproduction_size',              [],                @isnumeric); % 5
addOptional(p,'numberof_basisplanewaves',       -1,                @isnumeric); % 6
addOptional(p,'bright_weight',                  1.0,               @isnumeric); % 7
addOptional(p,'quiet_weight',                   2.5,               @isnumeric); % 8
addOptional(p,'unattended_weight',              0.05,              @isnumeric); % 9
addOptional(p,'brightzone_radius',              0.3,               @isnumeric); % 10
addOptional(p,'brightzone_pos_angle',           180,               @isnumeric); % 11
addOptional(p,'brightzone_pos_distance',        0.6,               @isnumeric); % 12
addOptional(p,'brightzone_geometry',            'circle',          @ischar);    % 13
addOptional(p,'brightzone_size',                [],                @isnumeric); % 14
addOptional(p,'brightzone_source_angle',        15,                @isnumeric); % 15
addOptional(p,'brightzone_source_dist',         1.0,               @isnumeric); % 16
addOptional(p,'brightzone_source_type',         'pw',              @ischar);    % 17
addOptional(p,'quietzone_radius',               0.3,               @isnumeric); % 18
addOptional(p,'quietzone_pos_angle',            0,                 @isnumeric); % 19
addOptional(p,'quietzone_pos_distance',         0.6,               @isnumeric); % 20
addOptional(p,'quietzone_geometry',             'circle',          @ischar);    % 21
addOptional(p,'quietzone_size',                 [],                @isnumeric); % 22
addOptional(p,'numberof_loudspeakers',          -1,                @isnumeric); % 23
addOptional(p,'loudspeaker_radius',             1.5,               @isnumeric); % 24
addOptional(p,'maximum_frequency',              8000,              @isnumeric); % 25
addOptional(p,'angleto_firstloudspeaker',       90,                @isnumeric); % 26
addOptional(p,'angleof_loudspeakerarc',         180,               @isnumeric); % 27
addOptional(p,'loudspeaker_model',              'Genelec 8010A',   @ischar);    % 28
addOptional(p,'loudspeaker_object',             [] );                           % 29
addOptional(p,'angleof_loudspeakerarrcentre',   180,               @isnumeric); % 30
addOptional(p,'loudspeaker_spacing',            0.01,              @isnumeric); % 31
addOptional(p,'speaker_array_type',             'circle',          @ischar);    % 32
addOptional(p,'room_size',                      [],                @isnumeric); % 33
addOptional(p,'dimensionality',                 2,                 @isnumeric); % 34
addOptional(p,'recordingtype',                  'simulated',                 @isnumeric); % 35
parse(p, settings{:});

setup = MultizoneSoundfieldSetup( ...       
    p.Results.resolution, ...                   1
    p.Results.frequency, ...                    2
    p.Results.reproduction_radius, ...          3
    p.Results.reproduction_geometry, ...        4
    p.Results.reproduction_size, ...            5
    p.Results.numberof_basisplanewaves, ...     6
    p.Results.bright_weight, ...                7
    p.Results.quiet_weight, ...                 8
    p.Results.unattended_weight, ...            9
    p.Results.brightzone_radius, ...            10
    p.Results.brightzone_pos_angle, ...         11
    p.Results.brightzone_pos_distance, ...      12
    p.Results.brightzone_geometry, ...          13
    p.Results.brightzone_size, ...              14
    p.Results.brightzone_source_angle, ...      15
    p.Results.brightzone_source_dist, ...       16
    p.Results.brightzone_source_type, ...       17
    p.Results.quietzone_radius, ...             18
    p.Results.quietzone_pos_angle, ...          19
    p.Results.quietzone_pos_distance, ...       20
    p.Results.quietzone_geometry, ...           21
    p.Results.quietzone_size, ...               22
    p.Results.numberof_loudspeakers, ...        23
    p.Results.loudspeaker_radius, ...           24
    p.Results.maximum_frequency, ...            25
    p.Results.angleto_firstloudspeaker, ...     26
    p.Results.angleof_loudspeakerarc, ...       27
    p.Results.loudspeaker_model, ...            28
    p.Results.loudspeaker_object, ...           29
    p.Results.angleof_loudspeakerarrcentre, ... 30
    p.Results.loudspeaker_spacing, ...          31
    p.Results.speaker_array_type, ...           32
    p.Results.room_size, ...                    33
    p.Results.dimensionality);%                 34

end

function setup = MultizoneSoundfieldSetup( ...
    res, ...     1
    f,...        2
    R,...        3
    reprGeom,... 4
    reprSize,... 5
    N,...        6
    Wb,...       7
    Wq,...       8
    Wu,...       9
    rb,...       10
    angb,...     11
    disb,...     12
    zoGeomb,...  13
    zoSizeb,...  14
    srcangb,...  15
    srcdistb,... 16
    srctypeb,... 17
    rq,...       18
    angq,...     19
    disq,...     20
    zoGeomq,...  21
    zoSizeq,...  22
    L,...        23
    Rl,...       24
    fmax,...     25
    phi,...      26
    phiL,...     27
    spkrmod,...  28
    spkrobj,...  29
    phiLcent,... 30
    spkrspace,...31
    arr_type,... 32
    roomSz,...   33
    dims)%       34
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
quiet.Dimensionality  = dims;
bright.Dimensionality = dims;
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
soundfield.Dimensionality = dims;
soundfield.Geometry  = reprGeom;
soundfield = soundfield.setReproRegionSize( reprSize );
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
setup.Dimensionality = dims;
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