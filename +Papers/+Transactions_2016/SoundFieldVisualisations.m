clc;
clear;
close all;
tic;

%%
% SYS = Current_Systems.loadCurrentSRsystem;
SYS = Current_Systems.IEEETransactions_System_G;

%%
setup = [SYS.Main_Setup(:);SYS.Masker_Setup(:)];
for s = 1:numel(setup)        
    setup(s).Multizone_Soundfield = setup(s).Multizone_Soundfield.createSoundfield('DEBUG');    
    setup(s) = setup(s).calc_Loudspeaker_Weights();    
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
    SYS.publication_info.axes_gap, ...
    SYS.publication_info.axes_margins_height, ...
    SYS.publication_info.axes_margins_width, ...
    'centimeters');
FontSize = 10;
FontName = 'Times';


for s = 1:numel(setup)
axes(ha(s));
ax=gca;
setup(s).plotSoundfield( mag2db(abs(F{s})), 'scientific_L12', realistic, details);

text(0,size(F{s},1),1e3,['(' char(64+s) ')'],...
    ... 'BackgroundColor',[1 1 1 0.7], ...
    'FontName',FontName,'FontSize',FontSize, ...
    'Interpreter','latex', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');
if s ~= 1
    ax.Title = [];
else
    ax.Title.String = ' ';
end

ax.XLabel.Interpreter = SYS.publication_info.Interpreter;
ax.YLabel.Interpreter = SYS.publication_info.Interpreter;

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
drawnow;
hCB = colorbar; 
hCB.Location = 'manual';
hCB.Units = 'points';
hCB.Label.Interpreter = SYS.publication_info.Interpreter;
hCB.TickLabelInterpreter = SYS.publication_info.Interpreter;
hCB.Ticks = interp1(1:length(caxis),caxis,linspace(1,length(caxis),6));
hCB.TickLabels = cellfun(@strrep, ...
    num2cell(num2str(linspace( ax.CLim(1), ax.CLim(2),6)' ),2), ...
    repmat({' '},6,1), repmat({''},6,1),'un',0);

hCB.Label.String = 'Magnitude (dB)';

ax.Units = 'centimeters'; hCB.Units = 'centimeters';

hCB.Position(1) = ...
    ax.Position(1) + ...
    ax.Position(3) + ...
    SYS.publication_info.axes_gap(2);
hCB.Position(2) = ...
    ax.Position(2);
hCB.Position(3) = ...
    SYS.publication_info.axes_gap(2);
hCB.Position(4) = ...
    ax.Position(4) + ...
    SYS.publication_info.axes_gap(1) + ...
    ax.Position(4); % Colorbar height


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



% C(end+1) = pow2db(setup(1).Acoustic_Contrast);
% E(end+1) = mag2db(setup(1).MSE_Bright);

%%
disp(['   Contrast: ' num2str(pow2db(setup(1).Acoustic_Contrast)) 'dB']);
disp(['        MSE: ' num2str(mag2db(setup(1).MSE_Bright)) 'dB']);
disp(['Attenuation: ' num2str(setup(1).Attenuation_dB(1)) 'dB (±' num2str(setup(1).Attenuation_dB(2)) 'dB)']);

%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script