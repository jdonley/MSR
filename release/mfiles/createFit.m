function [fitresult, gof] = createFit(NoiseLevel, Speech_Intelligibility, Type_of_Fit, Fit_Options)
%CREATEFIT(NOISELEVEL,SI_WC_BRIGHT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : NoiseLevel
%      Y Output: SI_WC_Bright
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( NoiseLevel, Speech_Intelligibility );

% Set up fittype and options.
if nargin < 3
    ft = fittype( 'smoothingspline' );
elseif nargin>=3
    ft = fittype(Type_of_Fit);
end

% Fit model to data.
if nargin <4
    [fitresult, gof] = fit( xData, yData, ft );
elseif nargin >=4
    [fitresult, gof] = fit( xData, yData, ft, Fit_Options );
end



