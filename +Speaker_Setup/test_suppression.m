%clc;
clear;
%close all;
tic;

%%
SYS = Current_Systems.loadCurrentSRsystem;

f = 1700;
c = 343;
f_ = logspace(log10(100),log10(8000),100);

hi_vals = [];talker=[];canceled=[];supp_confint=[];
%%
for fi = 1:numel(f_)
    f = f_(fi);
    setup = SYS.Main_Setup(1);
    %setup = SYS.Masker_Setup;
    
    setup.Multizone_Soundfield.Quiet_Zone = ...
        setup.Multizone_Soundfield.Quiet_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup.Multizone_Soundfield.Bright_Zone = ...
        setup.Multizone_Soundfield.Bright_Zone.setDesiredSoundfield(true, f, 'suppress_output');
    setup.Multizone_Soundfield = setup.Multizone_Soundfield.setN( -1 ); %Auto set
    setup.Multizone_Soundfield = setup.Multizone_Soundfield.createSoundfield('DEBUG');
    
    setup.Loudspeaker_Dimensions(2) = 1e-5;
    setup = setup.calc_Loudspeaker_Locations;
    setup = setup.calc_Loudspeaker_Weights();
    
    if strcmpi(setup.Speaker_Array_Type, '2line')
        [x,y]=pol2cart(setup.Loudspeaker_Locations(:,1),setup.Loudspeaker_Locations(:,2));
        d = sum(diff([x([1 end]),y([1 end])]).^2).^0.5;
        setup.Loudspeaker_Weights(1:end/2) = setup.Loudspeaker_Weights([end:-1:end/2+1])*exp(1i* (pi- d/(c/f) * 2*pi));
        %delSrcI = y > -0.1 & y < 0.1;
        %setup.Loudspeaker_Weights(delSrcI) = 0;
    end
    
    setup = setup.reproduceSoundfield('DEBUG');
    
    %setup.Loudspeaker_Locations(delSrcI,:) = nan;
    
    %%
    S_actual = setup.Bright_Samples;
    
    bz = setup.Multizone_Soundfield.Bright_Zone;
    mask = 1*bz.Soundfield_d_mask;
    mask(mask==0)=nan;
    S_desired = (bz.Soundfield_d.*mask);
    
    S_a_norm = S_actual ./ rms(abs(S_actual(~isnan(mask))));
    S_d_norm = S_desired ./ rms(abs(S_desired(~isnan(mask))));
    
    disp(fi);
    
    talker(fi,:) = abs( S_d_norm(~isnan(mask))  );
    canceled(fi,:) = abs( S_d_norm(~isnan(mask)) - S_a_norm(~isnan(mask)) );
    
end
%%
% hi=[];
% for i = 1:size(suppression,1)
%     figure(1);
%     h = histogram(mag2db(suppression(i,:)),[-45.05:0.1:0.05]);
%     hi(i,:) = h.Values;
% end
% figure(2)
% surf(f_,-45:0.1:0, mag2db(hi.' ./ size(suppression,2)),'linestyle','none');
% set(gca,'XScale','log');grid on;grid minor;
% xlim([90 9000]); ylim([-45 5]);view(2);

% dbConversion = 20/log(10);
% m = mean((suppression.'));
% v = var(suppression.');
% mu_ = log( m ./ sqrt( 1 + v./m.^2 ) ) *dbConversion;
% sigm = sqrt( log( 1 + v./m.^2 ) ) *dbConversion;


suppression = mag2db(canceled)-mag2db(talker);
mu_ = mean(suppression,2).';
sigm = sqrt( var(suppression,0,2) ).';


figure(3);
ha=area([f_; f_].'/1e3, [mu_+sigm; -2*sigm].','LineStyle','none','showbaseline','off'); hold on;
ha(1).FaceColor=[1 1 1]*1;drawnow;
ha(2).FaceColor=[1 1 1]*0;drawnow;
ha(1).FaceAlpha=0.0;
ha(2).FaceAlpha=0.5;
plot(f_/1e3,mu_,'color',[1 1 1]*0); hold on;

set(gca,'XScale','log');grid on;grid minor;
xlim([90 10000]/1e3); ylim([-45 5]);
xlabel('Frequency (kHz)');
ylabel('Attenuation (dB)');
hold off;
%%
%fprintf(Speaker_Setup.printSetupDetails(setup));

%%
tEnd = toc;
fprintf('Execution time: %dmin(s) %fsec(s)\n', floor(tEnd/60), rem(tEnd,60)); %Time taken to execute this script