%clc;
clear;
%close all;
tic;

%%
% SYS = Current_Systems.loadCurrentSRsystem;

SYS = Current_Systems.IEEETransactions_System_G;

%%
f = 1000;
%  f = Broadband_Tools.getAliasingFrequency(SYS.Main_Setup(2))*343/2/pi;
c = 343;
% Freqs = 200:100:3000;
C=[];E=[];
%  for f = Freqs
% fprintf('%.0f\n',f);
%%
setup = [SYS.Main_Setup(:);SYS.Masker_Setup(:)];

for s = 1:numel(setup)
    
%         setup(s).Multizone_Soundfield.Radius = 0.91;
%     setup(s).Multizone_Soundfield.UnattendedZ_Weight = 0;
    
%     setup(s).Multizone_Soundfield.Quiet_Zone = ...
%         setup(s).Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
%     setup(s).Multizone_Soundfield.Bright_Zone = ...
%         setup(s).Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
%     setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.setN( -1 ); %Auto set
%     setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createEmptySoundfield('DEBUG');
%     if setup(s).Loudspeaker_Count>1
        setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createSoundfield('DEBUG');
%     end
    
    setup(s) = setup(s).calc_Loudspeaker_Weights();
%     if s == 2
%        setup(s).Radius = (setup(s-1).Loudspeaker_Dimensions(1)*setup(s-1).Loudspeaker_Count/2)/2;
% 
%        setup(s).Multizone_Soundfield.Radius = setup(s).Radius;
%    end
    setup(s) = setup(s).reproduceSoundfield('DEBUG');
    
end


%%
if exist('fH'), if isvalid(fH), close(fH); end; end

figNums = [101,102,103];
realistic = false;
details.DrawDetails = false;
details.zoneLineWid = 1.0;
details.arrowLineWid = 0.5;
details.arrowLength = 3;
details.arrowAngle = 30;
details.arrowBuffer = 2;
details.lblFontSize = 10;

for s = 1:numel(setup)
 pk(s) = max(abs((setup(s).Bright_Samples(:))))*setup(s).res;
 
 F{s} = setup(s).Soundfield_reproduced*setup(s).res;
end


% clipFact = 2;
for s = 1:numel(setup)
gainNorm = 1/max(pk(s)); % Normalise to the maximum of all bright peaks
    F{s} = F{s}*gainNorm; pk(s) = pk(s)*gainNorm;
    
%     F{s}(abs(F{s})>clipFact*pk(s))=nan;
end


% close all;
fH = figure('Name',SYS.publication_info.FigureName);
dimSz = SYS.publication_info.subPlotDims;
ha = tightPlots( ...
    dimSz(2), ...
    dimSz(1), ...
    SYS.publication_info.figure_width, ...
    SYS.publication_info.axis_aspect_ratio, ...
    [0.5 -0.4], ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');
FontSize = 10;
FontName = 'Times';


for s = 1:numel(setup)
axes(ha(s));
ax=gca;
setup(s).plotSoundfield( mag2db(abs(F{s})), 'scientific_L12', realistic, details);

text(10,size(F{s},1)-FontSize-10,1e3,['(' char(64+s) ')'],...
    'BackgroundColor',[1 1 1 0.7], ...
    'FontName',FontName,'FontSize',FontSize)
ax.Title.String = '';

[c,r] = ind2sub(dimSz,s);
if c ~= 1
 ax.YLabel = [];
 ax.YTickLabel = [];
end
if r ~= dimSz(2)
 ax.XLabel = [];
 ax.XTickLabel = [];    
end

% clim_=[-1 1].*pk(s);
clim_=[-25 mag2db(abs(pk(s)))];
ax.CLim = clim_;
colorbar off;

end

% 
 hCB = colorbar; hCB.Units = 'points';
 hCB.Ticks = interp1(1:length(caxis),caxis,linspace(1,length(caxis),5));
 hCB.TickLabels = num2str(linspace( ax.CLim(1), ax.CLim(2),5)' );
 hCB.Label.String = 'Magnitude (dB)';

 
 tightfig;
 
 % 
% set(fH.Children, 'Units','Points')
% for c = 1:numel(fH.Children)
%  fH.Children(c).Position(2) = fH.Children(c).Position(2)+20;
% end



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



C(end+1) = pow2db(setup(1).Acoustic_Contrast);
E(end+1) = mag2db(setup(1).MSE_Bright);

%%
disp(['   Contrast: ' num2str(pow2db(setup(1).Acoustic_Contrast)) 'dB']);
disp(['        MSE: ' num2str(mag2db(setup(1).MSE_Bright)) 'dB']);
disp(['Attenuation: ' num2str(setup(1).Attenuation_dB(1)) 'dB (±' num2str(setup(1).Attenuation_dB(2)) 'dB)']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script