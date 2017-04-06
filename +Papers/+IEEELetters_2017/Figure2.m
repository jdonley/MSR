%clc;
%clear;
%close all;
tic;

%%
SYS = Current_Systems.IEEELetters2017_System_A;

%%
f = 1000;
c = SYS.signal_info.c;

%%
setup = [SYS.Main_Setup(:);SYS.Masker_Setup(:)];
for s = 1:numel(setup)
    
%         setup(s).Multizone_Soundfield.Radius = 0.91;
    setup(s).Multizone_Soundfield.UnattendedZ_Weight = 0;
    
    setup(s).Multizone_Soundfield.Quiet_Zone = ...
        setup(s).Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup(s).Multizone_Soundfield.Bright_Zone = ...
        setup(s).Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.setN( -1 ); %Auto set
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createSoundfield('DEBUG');
    
    setup(s) = setup(s).calc_Loudspeaker_Weights();
    if s == 2
%         setup(s).Radius = (setup(s-1).Loudspeaker_Dimensions(1)*setup(s-1).Loudspeaker_Count/2)/2;

%         setup(s).Multizone_Soundfield.Radius = setup(s).Radius;
   end
    setup(s) = setup(s).reproduceSoundfield('DEBUG');
    
end

%%
figNums = [101,102,103];
realistic = true;
details.DrawDetails = false;
details.zoneLineWid = 1.5;
details.arrowLineWid = 0.4;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 12;
details.NTicks = [5, 5]; % Number of ticks in X and Y
details.SYS = SYS;

pk(1) = max(abs(setup(1).Bright_Samples(:)));
pk(2) = max(abs(setup(2).Bright_Samples(:) + setup(3).Bright_Samples(:) ));


ZM = setup(1).Soundfield_reproduced*setup(1).res * SYS.Room_Setup.Wall_Reflect_Coeff;
ZT = setup(2).Soundfield_reproduced*setup(2).res;
ZI = setup(3).Soundfield_reproduced*setup(3).res * SYS.Room_Setup.Wall_Reflect_Coeff;

close all;

ha = tightPlots( 2, 1, ...
SYS.publication_info.figure_width, ...
SYS.publication_info.axis_aspect_ratio, ...
SYS.publication_info.axes_gap, ...
SYS.publication_info.axes_margins_height, ...
SYS.publication_info.axes_margins_width, ...
'centimeters');

FontSize = 9;
FontName = 'Times';
axes(ha(1));
ax(1)=gca;
setup(1).plotSoundfield( ZT + ZI, 'scientific_D1', realistic, details);
text(500-10,size(ZT,1)-FontSize*2-10,1e3,'(A)','FontName',FontName,'FontSize',FontSize)
ax(1).Title.String = '';%'Pressure Soundfield of Talker';
ax(1).XLabel = [];
ax(1).XTickLabel = [];
clim_=[-1 1].*pk(1)*setup(1).res;
ax(1).CLim = clim_;
colorbar off

axes(ha(2))
ax(2)=gca;
setup(1).plotSoundfield( (ZT + ZI) - ZM, 'scientific_D1', realistic, details);
text(500-10,size(ZT,1)-FontSize*2-10,1e3,'(B)','FontName',FontName,'FontSize',FontSize)
ax(2).Title=[];
ax(2).CLim=clim_;
colorbar off

for i = 1:2
    ax(i).FontSize=FontSize;
end
set(gcf,'color','w');

drawnow;
pause(1);
tightfig;

SYS.publication_info.FigureName = 'IEEE_Letters2017_1';
SYS.publication_info.print_fmt = 'png';
Publication.saveFigureForPublication( SYS, gcf );

%%
disp(['   Contrast: ' num2str(mag2db(setup(1).Acoustic_Contrast)) 'dB']);
disp(['        MSE: ' num2str(mag2db(setup(1).MSE_Bright)) 'dB']);
disp(['Attenuation: ' num2str(setup(1).Attenuation_dB(1)) 'dB (±' num2str(setup(1).Attenuation_dB(2)) 'dB)']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script