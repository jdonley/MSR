%clc;
clear;
%close all;
tic;

%%
SYS = Current_Systems.OculusPres_System_A;

ff = [...
    ...1000 ...
    1700 ...
    2500 ...
    5000];
c = 343;

setup = SYS.Main_Setup(ones(numel(ff),1));
for s = 1:numel(ff)
    f = ff(s);


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
try close('111'); catch; end

figNums = [101,102,103];
realistic = false;
details.DrawDetails = false;
details.zoneLineWid = 1.5;
details.arrowLineWid = 0.4;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 12;
details.NTicks = [5 5];

pk(1) = max(abs((setup(1).Bright_Samples(:))))*setup(1).res; % Masker (loudspeakers)
pk(2) = max(abs((setup(2).Bright_Samples(:))))*setup(2).res; % Talker
pk(3) = max(abs((setup(3).Bright_Samples(:))))*setup(3).res; % Image sources (reflections)

Z1 = setup(1).Soundfield_reproduced*setup(1).res;
Z2 = setup(2).Soundfield_reproduced*setup(2).res;
Z3 = setup(3).Soundfield_reproduced*setup(3).res;

gainNorm = 1/pk(1);
Z1 = Z1*gainNorm; pk(1) = pk(1)*gainNorm;
Z2 = Z2*gainNorm; pk(2) = pk(2)*gainNorm;
Z3 = Z3*gainNorm; pk(3) = pk(3)*gainNorm;

clipFact = 10;
Z1(abs(Z1)>clipFact*pk(1))=nan;
Z2(abs(Z2)>clipFact*pk(2))=nan;
Z3(abs(Z3)>clipFact*pk(3))=nan;



fH = figure(111);
ha = tightPlots( 1, 3, ...
    SYS.publication_info.figure_width*3, ...
    SYS.publication_info.axis_aspect_ratio, ...
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');

FontSize = 16;
FontName = 'Times';

%%%
axes(ha(1));
ax=gca;
setup(1).plotSoundfield( Z1, 'scientific_D1', realistic, details);
% text(10,size(Z1,1)-FontSize-10,1e3,'(A)',...
%     'BackgroundColor',[1 1 1 0.7],'FontName',FontName,'FontSize',FontSize)
ax.Title.String = '';%'Pressure Soundfield of Talker';
clim_=[-1 1].*pk(1);
ax.CLim = clim_;
colorbar off

%%%
axes(ha(2));
ax=gca;
setup(2).plotSoundfield( Z2, 'scientific_D1', realistic, details);
% text(10,size(Z2,1)-FontSize-10,1e3,'(B)',...
%     'BackgroundColor',[1 1 1 0.7],'FontName',FontName,'FontSize',FontSize)
ax.Title.String = '';%'Pressure Soundfield of Talker';
clim_=[-1 1].*pk(2);
ax.CLim = clim_;
colorbar off

%%%
axes(ha(3));
ax=gca;
setup(3).plotSoundfield( Z3, 'scientific_D1', realistic, details);
% text(10,size(Z3,1)-FontSize-10,1e3,'(C)',...
%     'BackgroundColor',[1 1 1 0.7],'FontName',FontName,'FontSize',FontSize)
ax.Title.String = '';%'Pressure Soundfield of Talker';
clim_=[-1 1].*pk(3);
ax.CLim = clim_;
colorbar off

drawnow();

tightfig;




%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script