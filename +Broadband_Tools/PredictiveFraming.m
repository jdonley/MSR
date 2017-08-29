function [ PredictedFrames, PredictedSignal ] = PredictiveFraming( Frames, Buffer, Npredict, method )
%This function accepts a framed signal and corresponding buffers to predict forward Npredict samples.
%
% Syntax:   [ PredictedFrames, PredictedSignal ] = PredictiveFraming( Frames, Buffer, Npredict, method )
%
% Inputs:
%	  Frames - The framed signal for which to predict ahead of each frame
%              by Npredict
%	  Buffer - The set of buffers corresponding to each frame of 'Frames'
%   Npredict - The number of samples to predict ahead of time at each frame
%     method - The autoregressive method to use for prediction
%
% Outputs:
% 	PredictedFrames - The set of frames containing the predicted samples
%   PredictedSignal - The frame based predicted signal
%
% Example:
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
%
% See also: predict, forecast, ar, lpc, aryule, arburg, arcov, armcov

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2016-2017
% Date: 16 August 2016
% Version: 0.1 (16 August 2016)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_p = Frames.'; % Frames_predicted
B = Buffer.';
[m,n1] = size(F_p);
[l,n2] = size(B);

if n1 ~= n2
    error('Number of frames mismatches with buffer.')
end

if strcmpi(method,'lpc') % Stable
    a = lpc(B,m-1);
elseif strcmpi(method,'yule') % Yule-Walker (yule) and Autocorrelation method of Least-Squares (lpc) are equivalent
    a = aryule(B,m-1);
elseif strcmpi(method,'burg') % Stable
    a = arburg(B,m-1);
elseif strcmpi(method,'cov') % Unstable
    a = arcov(B,m/2);
elseif strcmpi(method,'mcov') % Unstable
    a = armcov(B,m-1);
end

A = -[zeros(n1,1) a(:,2:end)].';
y = zeros(Npredict,n1);
for f=1:n1-1
    [~, zf] = filter(A(:,f), 1, B(:,f));
    y(:,f) = filter([0 0], -a(f,:), zeros(1, Npredict), zf);
end
y_p = [B((end-m+Npredict+1):end,:); ...
    y;] ; % Predicted part
y_p(isnan(y_p))=0;

F_p = [F_p(m/2+1:end,:); ... % Latest part-frame known values
    y_p(1:m/2,:); ...
    y_p];

PredictedFrames = reshape(F_p,m,[]).';
PredictedSignal = [zeros(m,1); y(:)];

%% Debug plots
% f=figure(12321);f.Name='PredictiveFraming'; hold off;
% for i = 0:size(Frames,1)
%     
%     subplot(3,1,1); xlim([-10*Npredict 0] + (i+1)*Npredict+m);  grid on;
%     plot((1:m)+i*Npredict,Frames(i+1,:),'linew',1.5);hold on;
%     
%     subplot(3,1,2); xlim([-10*Npredict 0] + (i+1)*Npredict+m);  grid on;
%     plot((-l+1:0)+m+i*Npredict,Buffer(i+1,:),'linew',1.5);hold on;
%     
%     subplot(3,1,3); xlim([-10*Npredict 0] + (i+1)*Npredict+m);  grid on;
%     plot((1:m)+Npredict+i*Npredict,PredictedFrames(2*i+2,:),'linew',1.5);hold on;
%     
%     pause(0.1);
%     drawnow;
%     0;
% %     input('enter to continue...');
% end



end

