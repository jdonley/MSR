clear; close all;
%%
SYS = Current_Systems.loadCurrentSRsystem;

f = 1000;

v = VideoWriter('Ripples.avi','Uncompressed AVI');
open(v);

%%
for SET = 1:numel(SYS.Main_Setup)
    setup = SYS.Main_Setup(SET);
    %setup = SYS.Masker_Setup;
    
    setup.Multizone_Soundfield.Quiet_Zone = ...
        setup.Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup.Multizone_Soundfield.Bright_Zone = ...
        setup.Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup.Multizone_Soundfield = setup.Multizone_Soundfield.setN( -1 ); %Auto set
    setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
    setup = setup.calc_Loudspeaker_Weights();
    setup = setup.reproduceSoundfield('DEBUG');
    
    %%
    figNums = [101,102,103];
    realistic = false;
    
    pk = max(abs((setup.Bright_Samples(:))));
    %pk = max(abs((setup.Quiet_Samples(:))));
    Z = setup.Soundfield_reproduced*setup.res;
    
    %figure(figNums(1))
    %setup.plotSoundfield( Z, 'scientific_D1A');
    
    
    %%
    view_pos = [40,40];
    Zscale = setup.res/10;
    depthMult = 3;
    depth = setup.res;
    field = real(Z);
    hf = figure(figNums(1));
    ha = surf(field/setup.res*Zscale);
    
    pk=max(abs(field(:)));
    fn= field/pk;
    % set(ha, ...
    %     'FaceAlpha','interp',...
    %     'AlphaData', (fn<0 -1) +1- fn.*(fn>=0) );
    axis equal
    ax = axis;
    axis([ax(1:4)+[-1 1 -1 1]*10 [-1 1]*depth] );
    ax = axis;
    caxis(ax(5:6));
    
    shading interp
    blue = linspace(0.4, 1.0, 25).' ; cm = [blue*0, blue*0, blue]; %'// create blue colormap
    colormap(cm);
    axis off
    lightangle(view_pos(1)-270,view_pos(2)) ;   %// add a light source
%     lightangle(view_pos(1)-180,view_pos(2)) ;   %// add a light source
    set(ha,'FaceLighting','phong','AmbientStrength',.3,'DiffuseStrength',.8,'SpecularStrength',.9,'SpecularExponent',25,'BackFaceLighting','unlit')
    view(view_pos(1),view_pos(2))
    set(gcf,'color',[1 1 1]*0.75);
    
    NsrcsTxt=['Sources: ' num2str(setup.Loudspeaker_Count)];
    DistTxt=['Distribution: ' upper(setup.Speaker_Array_Type(1)) lower(setup.Speaker_Array_Type(2:end))];
    % ylabel(NsrcsTxt);
    % title(NsrcsTxt);
    text(size(field,1)*0.6,size(field,2)*0.6,setup.res*1.8,NsrcsTxt,'FontName','arial','FontSize',20);
    text(size(field,1)*0.6,size(field,2)*0.6,setup.res*2.0,DistTxt,'FontName','arial','FontSize',20);
    
    hf.Units = 'normalized';
    hf.OuterPosition = [0 0 1 1];
%     tightfig(hf);
    
    %%
    spkrDeets.sphereRad = 5; %centimeters
    hSpkrs = setup.speakerPlots(field/setup.res*Zscale, false, spkrDeets);
    [sl_x, sl_y] = pol2cart(setup.Loudspeaker_Locations(:,1), (setup.Loudspeaker_Locations(:,2) * setup.res)-1);
    sl_x = ceil(sl_x + size(field,1)/2);
    sl_y = ceil(sl_y + size(field,2)/2);
    
    
    hold on;spkrZ=[];
    for i = 1:setup.Loudspeaker_Count
        spkrZ(i) = max(hSpkrs.sphere(i).ZData(:));
        hSpkrs.rod(i) = plot3([ax(1) sl_x(i)], ...
            [sl_y(i) sl_y(i)], ...
            [spkrZ(i) spkrZ(i)], ...
            'k','linewidth',2);
        hSpkrs.sphere(i).LineStyle = 'none';
        hSpkrs.sphere(i).FaceColor = 'r';
    end
    hSpkrs.rodLine = plot3(ones(1,numel(spkrZ))*ax(1), ...
            sl_y.', ...
            spkrZ, ...
            'k','linewidth',2);
    hold off;
    
    %%
    zones = [setup.Multizone_Soundfield.Bright_Zone, ...
            setup.Multizone_Soundfield.Quiet_Zone];
    for z_ = 1:numel(zones)
        zone_ = zones(z_);
        [cx,cy,cz] = cylinder(zone_.Radius_q,100);
        cx = ((cx + zone_.Origin_q.X)*setup.res + size(field,1)/2);
        cy = ((cy + zone_.Origin_q.Y)*setup.res + size(field,2)/2);
        cz = (~~cz-0.1)*depth*0.1;
        hold on;
        hcyl = surf( cx,cy,cz, ...
            'FaceColor',char('g'+(z_-1)*11), 'FaceAlpha','interp','AlphaData',flipud(~~(cz-min(cz(:)))), 'EdgeColor','none');
        hold off;
    end
    
    %%
    x = 1:size(field,1);
    y = 1:size(field,2);
    pcol = [0.2857,1,1.0] ;
    
    minLevel = -depth;
    xface = [ x x(end) x(1) ] ;
    yface = [ y y(end) y(1) ] ;
    face0 = zeros(size(xface)) ;
    faceZ = [zeros(size(x)) minLevel minLevel] ;
    BptsX = [x([2 end-1]) y([end-1 2])];
    BptsY = circshift(BptsX,[0,1]);
    faceB = ones(1,numel(BptsX)) * minLevel;
    
    pt(1) = handle( patch( xface , face0+1             , faceZ ,pcol ) ) ;
    pt(3) = handle( patch( xface , face0+size(field,2) , faceZ ,pcol ) ) ;
    pt(2) = handle( patch( face0+1             , yface , faceZ ,pcol ) ) ;
    pt(4) = handle( patch( face0+size(field,1) , yface , faceZ ,pcol ) ) ;
    ptB   = handle( patch( BptsX , BptsY       , faceB ,pcol ) ) ;
    
    pt = handle(pt) ;
    ptB = handle(ptB) ;
    set(pt,  'Facecolor',pcol , 'FaceAlpha',0.5 , 'EdgeColor','none')
    set(ptB, 'Facecolor',pcol*0.5 , 'FaceAlpha',1.0 , 'EdgeColor','none')
    
    
    %%
    for a = 0:20:360*10
        sf = real(Z.*exp(-1i*a/180*pi))/setup.res*Zscale; %shifted field -> sf
        
        set( ha, 'ZData', sf );
        
        %     hFtmp=figure;
        %     ha_alph = surfl(1:size(sf,1),1:size(sf,2),sf/setup.res*Zscale,view_pos+[-180 0]);
        %     CData_ = (ha_alph.CData - min(ha_alph.CData(:)));
        %     CData_ = 0.5 + CData_ / max(CData_(:))*2*0.5;
        %     CData_(CData_>1)=1;
        %     CData_(1)=1;
        %     CData_(end)=0;
        %     set(ha, ...
        %     'FaceAlpha','interp',...
        %     'AlphaData', CData_ );
        %     hFtmp.delete;
        
        Z_ = get( ha,'ZData') ;
        Z_(isnan(Z_))=0;
        pt(1).ZData(1:end-2) = Z_(1,:) ;
        pt(2).ZData(1:end-2) = Z_(:,1) ;
        pt(3).ZData(1:end-2) = Z_(end,:) ;
        pt(4).ZData(1:end-2) = Z_(:,end) ;
        
        spkrZvals = sf(sub2ind(size(sf),sl_y,sl_x));
        for i = 1:setup.Loudspeaker_Count
            hSpkrs.sphere(i).ZData =  hSpkrs.sphere(i).ZData-min(hSpkrs.sphere(i).ZData(:))+spkrZvals(i) ;            
            spkrZ(i) = max(hSpkrs.sphere(i).ZData(:));
            hSpkrs.rod(i).ZData = repmat(spkrZ(i),1,2);
        end
        hSpkrs.rodLine.ZData = spkrZ;
        
        drawnow;
%         pause(0.001);
        set(gcf,'Renderer','zbuffer');
        writeVideo(v,getframe);
    end
    
end

%%
close(v);