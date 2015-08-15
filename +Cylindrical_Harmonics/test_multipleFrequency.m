clc;
clear;
close all;

%%
sp1 = spatial_zone(500, 0.5, 32, 'pw');

sp2 = spatial_zone(1600, 0.5, -40, 'pw');

%%
gs1 = global_soundfield;
gs1 = gs1.addSpatialZone(sp1, 0.5, 32);

gs2 = global_soundfield;
gs2 = gs2.addSpatialZone(sp2, 2.0, -40);

%%
gs1 = gs1.createGlobalSoundfield(gs2.getRadius_p_FromZones);

gs2 = gs2.createGlobalSoundfield();


%%
gs2.GlobalSoundfield_d = gs2.GlobalSoundfield_d + gs1.GlobalSoundfield_d;
gs2.plotSoundfield();