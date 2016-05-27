function [ Hrz_Vec, Res_Matrix, Res_trend, Res_area, Res_CI, CI_vec ] = generatePlotData( X_vals,Res,ConfInt_Low,ConfInt_Up, Type_of_Fit, X_vals_buffer )
%GENERATEPLOTDATA Organises input data into nice trending plottable data
%for plots, area plots and errorbar plots.
% 
% Syntax:	[Hrz_Vec, Res_Matrix, Res_trend, Res_area, Res_CI, CI_vec] = GENERATEPLOTDATA(X_vals,Res,ConfInt_Low,ConfInt_Up) Explain usage here
% 
% Inputs: 
% 	input1 - Description
% 	input2 - Description
% 	input3 - Description
% 
% Outputs: 
% 	output1 - Description
% 	output2 - Description
% 
% Example: 
% 	Line 1 of example
% 	Line 2 of example
% 	Line 3 of example
% 
% See also: List related files here

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2015
% Date: 02 September 2015 
% Revision: 0.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    Type_of_Fit =  'smoothingspline' ;
end
if nargin < 6
   X_vals_buffer = [0 0]; 
end
    X_absdiff = abs(diff(X_vals));
    X_absdiff(X_absdiff == 0)=[];
    Num_Audio_Files = max(histcounts(X_vals,-100:min(X_absdiff):100));
            
% Order Results
    [~,I] = sort(X_vals);
    X_vals=X_vals(I);
    Res=Res(I);
    ConfInt_Low=ConfInt_Low(I);
    ConfInt_Up=ConfInt_Up(I);
    
    
% Create Matrices of Results
    %Masking noise Level
    NoiseLevel_Matrix = reshape(X_vals,Num_Audio_Files,size(X_vals,1)/Num_Audio_Files);
    Hrz_Vec = double(NoiseLevel_Matrix(1,:));
    
    %Speech Intelligibility Results
    Res_Matrix = reshape(Res,Num_Audio_Files,size(Res,1)/Num_Audio_Files);
    
    % Confidence Intervals from Spatial Sampling Points
    ConfInt_Low_M = reshape(ConfInt_Low,Num_Audio_Files,size(ConfInt_Low,1)/Num_Audio_Files);
    ConfInt_Up_M  = reshape(ConfInt_Up ,Num_Audio_Files,size(ConfInt_Up ,1)/Num_Audio_Files);
        
    %Average Spatial Sampling Confidence Intervals
    ConfInt_M = [mean(ConfInt_Low_M,1)' mean(ConfInt_Up_M,1)'];
    
        %% Calculate confidence intervals
    Res_CI = Tools.confidence_intervals(Res_Matrix, 95);
    
    %% Fit a trendline to the results
    Fit_Options = fitoptions( 'Method', Type_of_Fit );
    if strcmp(Type_of_Fit,'smoothingspline')
        Fit_Options.SmoothingParam = 1.0;
    end
    [Res_trend, ~] = Results.createFit(double(X_vals),Res_Matrix(:), Type_of_Fit, Fit_Options);
    
    %Trendline plot vectors
    CI_vec= linspace(Hrz_Vec(1)-X_vals_buffer(1),Hrz_Vec(end)+X_vals_buffer(2),100);
    
    temp1 = repmat(Res_CI(:,1)',Num_Audio_Files,1); temp1 = temp1(:);
    temp2 = repmat(Res_CI(:,2)',Num_Audio_Files,1); temp2 = temp2(:);
    if all((numel(temp1) ~= 1) & (numel(temp2) ~= 1))
        CI_trend = struct('Upper', ...
            Results.createFit(double(X_vals),temp1,Type_of_Fit,Fit_Options), ...
            'Lower', ...
            Results.createFit(double(X_vals),temp2,Type_of_Fit,Fit_Options));        
        
        % Calculate areas
        Res_area = ([(Res_trend(CI_vec)) + CI_trend.Upper(CI_vec), ...
            CI_trend.Lower(CI_vec) - CI_trend.Upper(CI_vec)]);
        
    else
        CI_trend = [];
        Res_area = [];
    end
    


end

