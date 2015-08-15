%% Plot Time-Frequency Component idea.
Input_file_path = '+Miscellaneous\+Speech_Files\';
Input_file_ext  = '.wav';
files = dir([Input_file_path '*' Input_file_ext]);
file = files';
file = file(1);

[Z, Frequencies_] = Broadband_Tools.FFT_custom([Input_file_path file.name(1:end-4) Input_file_ext], 1024, 16000, 0.5);
Z = Z(1:1:end,1:1:end);
Z = Z ./ max(abs(Z(:)));
f=Frequencies_(1:1:end);

Z=abs(Z);

% h=bar4viacolor(struct('Z',Z));
% tick = get(gca,'XTick');
% set(h,'edgecolor','none');
% set(gca,'XTickLabel',f(tick));
% %set(gca,'XScale','log');
% caxis([0 0.1]);

h=surf(f, 1:size(Z,1), Z, 'LineStyle','none');
set(gca,'XScale','log');
view(2);
caxis([0 0.1]);