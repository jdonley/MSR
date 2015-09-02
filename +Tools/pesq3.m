function [ res ] = pesq3( reference_sig, degraded_sig, Fs )

temp_path = [pwd filesep '+Miscellaneous\+Temporary\'];
ref_path = [ 'tmp_ref.wav'];
deg_path = [ 'tmp_deg.wav'];

if ~exist(temp_path,'dir'); mkdir(temp_path); end

max_val = max(abs([reference_sig(:); degraded_sig(:)]));

audiowrite([temp_path ref_path], reference_sig / max_val, Fs);
audiowrite([temp_path deg_path], degraded_sig / max_val, Fs);

res = pesq2_mtlb(ref_path, ...
                 deg_path, ...
                 Fs, 'wb', [pwd filesep '+Tools\pesq2.exe'], ...
                 temp_path);

%delete([temp_path ref_path],[temp_path deg_path]);

end

