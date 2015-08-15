%Minimum frequency for accurate control
fmin = 1.25;
%Maximum frequency for accurate control
fmax=8;
%Maximum level of attenuation conrol
Lmax = 50;

%% Threshold in Quiet
frequencies = 0.01:0.01:20;
thresh = Perceptual_Tools.Threshold_in_Quiet(frequencies);


hold on;
%%
y = [0 0];
plot ([0.01 20],y,'--k','LineWidth',3);
set(gca,'XScale','log'); grid on;
ylabel('dB'); xlabel('Frequency');
axis([0.01 20 -20 120]);

x = [fmin fmin];
plot (x,[0 Lmax],'black','LineWidth',3);

x = [fmax fmax];
plot (x,[0 Lmax],'black', 'LineWidth',3);

y = [Lmax Lmax];
plot ([fmin fmax],y,'black', 'LineWidth',3);


% Masker
fm=-Inf;%kHz
Lm = 90;%dB
x = [fm fm];
plot (x,[0 Lm],'red', 'LineWidth',10);

%% The Bark Scale
f_maskee = 10:10:20000;
dz = Perceptual_Tools.FrequencyToBark( f_maskee ) - Perceptual_Tools.FrequencyToBark( fm*1000 );


%% Masking Curves and Spreading Functions
z_maskee = Perceptual_Tools.FrequencyToBark( f_maskee );
z_masker = Perceptual_Tools.FrequencyToBark( fm*1000  );

%ISO/IEC MPEG Psychoacoustic Model 2 speading function
mask = Perceptual_Tools.Masking_Curves.ISO_IEC_MPEG_Model2( z_maskee, z_masker, Lm);
Mask_Final = Perceptual_Tools.Combine_ThreshMask(thresh, mask);

%Modified Schroeder Spreading Function







%%
plot( frequencies, Mask_Final, 'LineWidth', 3);
Mask_Final(Mask_Final > Lmax) = Lmax;
Mask_Final(frequencies > fmax) = 0;
Mask_Final(frequencies < fmin) = 0;
H = area( frequencies, Mask_Final);
h=get(H,'children');
set(h,'FaceAlpha',0.5,'FaceColor',[0.5 1 0.5]);

%%
hold off;


