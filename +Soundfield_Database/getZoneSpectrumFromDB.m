function [ spectrum, frqs, DBfrqs ] = getZoneSpectrumFromDB( Zone, setup, LUT_resolution, Drive, fs )
%GETZONESPECTRUM Summary of this function goes here
%   Detailed explanation goes here

%% First, Load the relevant look-up tables
    method = {'new4', 'new3', 'new2', 'new'};
    for m = 1:length(method)
            method_ = [method, {'old_zones_swapped'}];
            [DB,err] = Soundfield_Database.loadDatabaseFromSetup( setup, LUT_resolution, Drive, method_{m} );
        if ~err
            break;
        end
    end

    Frequencies = DB.Frequencies;
    Weights = DB.Weights;
    
    %% Find correct frequencies for given sampling frequency    
    freqs1 = linspace(0, fs/2, fs/2 + 1);
    freqs1_ = freqs1(freqs1>=min(Frequencies) & freqs1<=max(Frequencies));    
    
    %% Ideal quiet zone levels are found with maximum contrast
    LUT_MagDiff = DB.Acoustic_Contrast__Weight_Vs_Frequency;%DB.Bright_Sample__Weight_Vs_Frequency - DB.Quiet_Sample__Weight_Vs_Frequency;
    
    % Noise weights
    LUT_MagDiff_interp = interp2(Frequencies,Weights,LUT_MagDiff,freqs1_.',Weights,'spline');
%     [~,I]=max(LUT_MagDiff_interp);
    weights = repmat(Weights(end),1,numel(freqs1_));
    
    %% Find spectrum from given weights
    if strcmpi(Zone,'Bright')
       DBSpects = DB.Bright_Sample__Weight_Vs_Frequency; 
    elseif strcmpi(Zone, 'Quiet')
       DBSpects = DB.Quiet_Sample__Weight_Vs_Frequency;
    end
	spectrum = Tools.interpVal_2D(DBSpects, Frequencies, Weights, freqs1_, weights, 'spline');    

    %% Ramp spectrum in and out if needed
    spectrum = [ ...
        linspace(            1,spectrum(1), length(freqs1(freqs1<min(Frequencies))) ), ...
        spectrum(:)', ...
        linspace(spectrum(end),          1, length(freqs1(freqs1>max(Frequencies))) )];

    %% logarithmically spaced samples
%     frqs = [0, logspace(log10(freqs1(2)), log10(fs/2), fs/2 )];
%     spectrum = interp1( freqs1, spectrum, frqs );
      frqs = freqs1;
      
      DBfrqs = Frequencies;
end

