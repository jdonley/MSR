bright = Orthogonal_Basis_Expansion.spatial_zone(2000, 0, 0.30, 'pw', 1.0, 0);
bright.res = 2000;  
bright = bright.setDesiredSoundfield(true);

m = bright.Soundfield_d_mask*1;
m(m==0) = nan;
bright.Soundfield_d = bright.Soundfield_d.*m * bright.res * 0.3*2 * 0.1;
%bright.Soundfield_d(imag(bright.Soundfield_d) > 0.0) = nan;
h = figure(1);
ax = surf(imag(bright.Soundfield_d),'LineStyle','none');
drawnow; pause(0.05);
colormap([0.8 0.8 0.8; ...
          1 1 1;]);      
colormap([0.0 0.0 0.0; ...
          0.2 0.2 0.2;]);
axis off equal;
ax(1).FaceAlpha = 1.0;
h.Color=[0 0 0];
h.Color=[1 1 1];
view([ 180 30 ]);