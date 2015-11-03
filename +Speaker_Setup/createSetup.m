function [ setup ] = createSetup( settings )
%CREATESETUP Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
p.CaseSensitive = false;

addOptional(p,'resolution',                     100,               @isnumeric);
addOptional(p,'frequency',                      2000,              @isnumeric);
addOptional(p,'reproduction_radius',            1.0,               @isnumeric);
addOptional(p,'numberof_basisplanewaves',       -1,                @isnumeric);
addOptional(p,'bright_weight',                  1.0,               @isnumeric);
addOptional(p,'quiet_weight',                   2.5,               @isnumeric);
addOptional(p,'unattended_weight',              0.05,              @isnumeric);
addOptional(p,'brightzone_radius',              0.3,               @isnumeric);
addOptional(p,'brightzone_pos_angle',           180,               @isnumeric);
addOptional(p,'brightzone_pos_distance',        0.6,               @isnumeric);
addOptional(p,'brightzone_source_angle',        15,                @isnumeric);
addOptional(p,'quietzone_radius',               0.3,               @isnumeric);
addOptional(p,'quietzone_pos_angle',            0,                 @isnumeric);
addOptional(p,'quietzone_pos_distance',         0.6,               @isnumeric);
addOptional(p,'numberof_loudspeakers',          -1,                @isnumeric);
addOptional(p,'loudspeaker_radius',             1.5,               @isnumeric);
addOptional(p,'maximum_frequency',              8000,              @isnumeric);
addOptional(p,'angleto_firstloudspeaker',       90,                @isnumeric);
addOptional(p,'angleof_loudspeakerarc',         180,               @isnumeric);
addOptional(p,'loudspeaker_model',              'Genelec 8010A',   @ischar);
addOptional(p,'angleof_loudspeakerarrcentre',   180,               @isnumeric);
addOptional(p,'loudspeaker_spacing',            0.01,              @isnumeric);
addOptional(p,'speaker_array_type',             'circle',          @ischar);

parse(p, settings{:});

setup = MultizoneSoundfieldSetup( ...
    p.Results.resolution, ...
    p.Results.frequency, ...
    p.Results.reproduction_radius, ...
    p.Results.numberof_basisplanewaves, ...
    p.Results.bright_weight, ...
    p.Results.quiet_weight, ...
    p.Results.unattended_weight, ...
    p.Results.brightzone_radius, ...
    p.Results.brightzone_pos_angle, ...
    p.Results.brightzone_pos_distance, ...
    p.Results.brightzone_source_angle, ...
    p.Results.quietzone_radius, ...
    p.Results.quietzone_pos_angle, ...
    p.Results.quietzone_pos_distance, ...
    p.Results.numberof_loudspeakers, ...
    p.Results.loudspeaker_radius, ...
    p.Results.maximum_frequency, ...
    p.Results.angleto_firstloudspeaker, ...
    p.Results.angleof_loudspeakerarc, ...
    p.Results.loudspeaker_model, ...
    p.Results.angleof_loudspeakerarrcentre, ...
    p.Results.loudspeaker_spacing, ...
    p.Results.speaker_array_type);

end

function setup = MultizoneSoundfieldSetup(res, f, R, N, Wb, Wq, Wu, rb, angb, disb, PWangb, rq, angq, disq, L, Rl, fmax, phi, phiL, spkrmod, phiLcent, spkrspace, arr_type)
%%
c = 343;
k = (f/c)*2*pi;
M = ceil(k*R);
kmax = (fmax/c)*2*pi;
Mmax = ceil(kmax*R);

%%
quiet  = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rq, 'quiet');
bright = Orthogonal_Basis_Expansion.spatial_zone(f, 0, rb, 'pw', 1.0, PWangb);
quiet.res  = res;
bright.res = quiet.res;
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
setup = setup.addMultizone_Soundfield(soundfield);
setup = setup.setRadius( Rl );

setup.Speaker_Arc_Angle = phiL;
setup.Angle_FirstSpeaker = phi;

if L < 1
    L = ceil( setup.Speaker_Arc_Angle/360 * 2*Mmax + 1 );
end
setup.Loudspeaker_Count = floor(L);

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