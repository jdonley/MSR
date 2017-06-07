%clc;
clear;
%close all;
tic;

%%
SYS = Current_Systems.loadCurrentSRsystem;
% SYS = Current_Systems.ICASSP2017_System_A;

f = 1000;
c = 343;

C=[];E=[];
% for f = 200:100:2000

%%
setup = [SYS.Main_Setup(:);SYS.Masker_Setup(:)];
for s = 1:3
    
%         setup(s).Multizone_Soundfield.Radius = 0.91;
%     setup(s).Multizone_Soundfield.UnattendedZ_Weight = 0;
    
    setup(s).Multizone_Soundfield.Quiet_Zone = ...
        setup(s).Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup(s).Multizone_Soundfield.Bright_Zone = ...
        setup(s).Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.setN( -1 ); %Auto set
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createEmptySoundfield('DEBUG');
    if setup(s).Loudspeaker_Count>1
        setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createSoundfield('DEBUG');
    end
    
    setup(s) = setup(s).calc_Loudspeaker_Weights();
%     if s == 2
%        setup(s).Radius = (setup(s-1).Loudspeaker_Dimensions(1)*setup(s-1).Loudspeaker_Count/2)/2;
% 
%        setup(s).Multizone_Soundfield.Radius = setup(s).Radius;
%    end
    setup(s) = setup(s).reproduceSoundfield('DEBUG');
    
end

%%
figNums = [101,102,103];
realistic = false;
details.DrawDetails = false;
details.zoneLineWid = 1.5;
details.arrowLineWid = 0.4;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 12;

 pk(1) = max(abs((setup(1).Bright_Samples(:))))*setup(1).res; % Masker (loudspeakers)
 pk(2) = max(abs((setup(2).Bright_Samples(:))))*setup(2).res; % Talker
 pk(3) = max(abs((setup(3).Bright_Samples(:))))*setup(3).res; % Image sources (reflections)
 
ZM = setup(1).Soundfield_reproduced*setup(1).res;
ZT = setup(2).Soundfield_reproduced*setup(2).res;
ZI = setup(3).Soundfield_reproduced*setup(3).res;

ZM(abs(ZM)>3*pk(1))=nan;
ZT(abs(ZT)>3*pk(2))=nan;
ZI(abs(ZI)>3*pk(3))=nan;
% Z2 = angle(Z);
% Z3 = abs(Z/setup.res);
% Z_ = mag2db((Z)./pk);

% close all;
figure(111);
ha = tightPlots( 3, 1, ...
SYS.publication_info.figure_width*2, ...
SYS.publication_info.axis_aspect_ratio, ...
SYS.publication_info.axes_gap, ...
SYS.publication_info.axes_margins_height, ...
SYS.publication_info.axes_margins_width, ...
'centimeters');

FontSize = 16;
FontName = 'Times';
axes(ha(1));
ax=gca;
setup(1).plotSoundfield( ZI, 'scientific_D1', realistic, details);
text(10,size(ZT,2)-FontSize/2-10,1e3,'(A)','FontName',FontName,'FontSize',FontSize)
ax.Title.String = '';%'Pressure Soundfield of Talker';
ax.XLabel = [];
ax.XTickLabel = [];
clim_=[-1 1].*pk(3);
ax.CLim = clim_;
colorbar off

axes(ha(2))
ax=gca;
setup(2).plotSoundfield( -ZM, 'scientific_D1', realistic, details);
text(10,size(ZT,2)-FontSize/2-10,1e3,'(B)','FontName',FontName,'FontSize',FontSize)
ax.Title=[];
ax.CLim=clim_;
colorbar off

axes(ha(3))
ax=gca;
setup(2).plotSoundfield( ZI-ZM, 'scientific_D1', realistic, details);
text(10,size(ZT,2)-FontSize/2-10,1e3,'(B)','FontName',FontName,'FontSize',FontSize)
ax.Title=[];
ax.CLim=clim_;
colorbar off

tightfig;


% figure(figNums(2)); hold off
% setup.plotSoundfield( Z2, 'scientific_L9', realistic, details);
% figure(figNums(3)); hold off
% setup.plotSoundfield( Z3, 'scientific_L9', realistic, details);



% for fn = 1:numel(figNums)
%     figure(figNums(fn));
%     if setup.Loudspeaker_Count > 1
%         R = [0 1].*size(Z,1) ; xlim(R);ylim(R);
%     else
%         R = [0 1].*size(Z,1) - size(Z,1) - x_*setup.res; xlim(R);ylim(R);
%     end
% end
%  caxis([-30, 0] );


%title('Small Zone Weight');

%  hold on;

%  Z = setup.Soundfield_reproduced;
%  %Z_ = -QualityGuidedUnwrap2D_r1(Z);
%  Z_ = GoldsteinUnwrap2D_r1(Z);
%
%  Zs = abs(Z);
%  Zs = Zs(1:10:end,1:10:end);
%  Z__ = Z_(1:10:end,1:10:end);
%
%  [U,V] = gradient( Z__(2:end-1,2:end-1) );
%  [X,Y] = meshgrid( 1:size(Z,1) , 1:size(Z,2) );
%  X = X(1:10:end,1:10:end);
%  Y = Y(1:10:end,1:10:end);
%  X = X(2:end-1,2:end-1);
%  Y = Y(2:end-1,2:end-1);
%  quiver3( X , Y, ones(size(U))*4, U .* abs(Zs(2:end-1,2:end-1)) , V .* abs(Zs(2:end-1,2:end-1)), zeros(size(U)), 1, 'k' );
%quiver3( X , Y, ones(size(U))*4, U  , V , zeros(size(U)), 1, 'k' );
%

C(end+1) = mag2db(setup(1).Acoustic_Contrast);
E(end+1) = mag2db(setup(1).MSE_Bright);

% end
%%
% figure(1010)
% hold on;
% plot(200:100:2000,C); hold off
% figure(1011)
% hold on;
% plot(200:100:2000,E); hold off
%%
disp(['   Contrast: ' num2str(mag2db(setup(1).Acoustic_Contrast)) 'dB']);
disp(['        MSE: ' num2str(mag2db(setup(1).MSE_Bright)) 'dB']);
disp(['Attenuation: ' num2str(setup(1).Attenuation_dB(1)) 'dB (�' num2str(setup(1).Attenuation_dB(2)) 'dB)']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script