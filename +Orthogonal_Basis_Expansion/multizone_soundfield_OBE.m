classdef multizone_soundfield_OBE
    % MULTIZONE_SOUNDFIELD_OBE - This class represents a multi-zone soundfield
    %                     using orthogonal basis expansion for reproduction 
    %                     by any single zone reproduction technique.
    %
    %   This class was implemented based on the publication:
    %   JIN, W., KLEIJN, W. B. & VIRETTE, D. Multizone soundfield reproduction 
    %   using orthogonal basis expansion.  Acoustics, Speech and Signal Processing 
    %   (ICASSP), 2013 IEEE International Conference on, 2013. IEEE, 311-315.
    %
    %
    %   Author: Jacob Donley, University of Wollongong, Australia
    %   Email: Jacob.Donley089@uowmail.edu.au
    %
    
    properties        
        % Settings
        res = 50;  % samples per metre %Resolution of soundfield
        Dimensionality = 2;          % 2D or 3D
        BrightZ_Weight      = 1.0;   % Bright Zone relative importance weight
        QuietZ_Weight       = 2.5;   % Quiet Zone relative importance weight
        UnattendedZ_Weight  = 0.05;  % Unattended Zone relative importance weight        
        Radius = 1.0;
        k_global = 2000 /343*2*pi;  % Frequency in wavenumber
        ModeLimitFactor = 1.0;      % Increase/Decrease the Mode Limit by a certain factor (default results in 16% error)
        Quiet_Zone;                 % spatial_zone object for the quiet zone
        Bright_Zone;                % spatial_zone object for the bright zone
        Geometry = 'circle';        % Geometry of the reproduction region
        ReproRegionSize = [];
        
        % Results
            % 3D
        Soundfield_desired = [];    % The complex values of the global desired sound field that is produced from this class.
        Soundfield_desired_mask = [];
        Quiet_Field = [];           % The complex values of the resultant quiet field only (normalised to the bright field)
        Bright_Field = [];          % The complex values of the resultant bright field only (normalised)        
        Bright_Error_Field = [];    % The error in the bright field
        Quiet_Error_Field = [];
            % 2D
        Alpha_Coeffs = [];           % Alpha_Coeffs are the global soundfield coefficients   
        Sound_sample_Avg = [];  %The average sound wave across the bright field in the direction of the desired planewave propagation
        Binaural_sample = [];   %TODO: Create binaural sample 
            % 1D
        Angle_Bright_Field = 0;
        Angle_Quiet_Field = 0;
        MSE_Bright_Field = 0;
        Err_dB_Bright_Field = 0;
        MSE_Quiet_Field = 0;
        Err_dB_Quiet_Field = 0;
        Acoustic_Brightness_Contrast = 0;
        Acoustic_Brightness_Contrast_dB = 0;
        R = [];
    end
    
    properties (Access = public)        
        N = 80;
        Delta_phi;
        F_angles = [];
        F = [];
        G = [];
        C = [];
        P = []; 
        w = [];
        Sd= [];
    end
    
%% Contructor
    methods 
        function obj = multizone_soundfield_OBE(~)
            obj.Delta_phi = 2 * pi / obj.N;
        end
    end

%% Private Methods
    methods (Access = public) %(Access = private) % usually private
            
        function obj = createEmptySoundfield(obj, Debug)
            if nargin < 2
                Debug = '';
            end
            %             if strcmp(Debug, 'Bright')
            %                 width = int16( ceil(obj.N^(1/2)) + ~mod(ceil(obj.N^(1/2)),2) );
            %                 obj.Soundfield_desired = zeros(width,width);
            %             else
            if strcmpi( obj.Geometry, 'circle' )
                width = int16(obj.res * obj.Radius * 2);
                height = width;
            elseif contains( lower(obj.Geometry), 'rect' )
                width  = obj.ReproRegionSize(1) * obj.res;
                height = obj.ReproRegionSize(2) * obj.res;
            end
            obj.Soundfield_desired = zeros(width,height);
            obj = obj.calcDesiredMask();
            %             end
        end
        
    end
    
%% Public Methods
    methods
        
%%%%%%%%%%%%%%%%%%  Main Soundfield Creation Function  %%%%%%%%%%%%%%%%%%%%

        function obj = createSoundfield(obj, Debug, Radius)
            if nargin < 2
                Debug = '';
            elseif nargin >= 3
               obj.Radius = Radius; 
            end            
            obj = obj.setWavenumberFromChildZone(); %Incase the child spatial zones have been changed
            
            M0 = obj.getGlobalModeLimit; 
            
            obj = obj.createEmptySoundfield(Debug);
            
            obj = obj.calc_Alpha_Coeffs(Debug);
                        
            if ~contains(Debug,'DEBUG') %% && ~strcmp(Debug, 'Bright') && ~strcmp(Debug, 'Quiet')
                obj = obj.createEmptySoundfield;
                O = length(obj.Soundfield_desired) / 2; %Set the centre of the zone for indexing
                theta = 0;
                r = 0;
                J = 0;
                e = 0;        
                M = -M0:M0;
                
                n=0;
                if ~contains(Debug,'suppress_output')
                    fprintf('--- Creating Global Soundfield ---\n');
                    fprintf('\tQuiet and Bright zones at %.0fHz.\n\tCompletion: ', obj.getFrequency);
                end
                
                for x = 1:O*2
                    for y = 1:O*2
                        theta = atan2(O - y, O - x);
                        r = ((O - x)^2 + (O - y)^2)^(1/2);
                        
                        r = r / obj.res;  %Change to spatial radius (/metre)
                        J = besselj(M, obj.k_global * r );
                        e = exp(1i * M * (-theta+pi/2));
                        obj.Soundfield_desired(x, y) = obj.Soundfield_desired(x, y) + sum(obj.Alpha_Coeffs .* J .* e); % Sum pressures
                        
                        
                    end
                    if ~contains(Debug,'suppress_output')
                        fprintf(repmat('\b',1,n));
                        n=fprintf('%.2f%%', x / (O*2) * 100);
                    end
                end
                
                if ~contains(Debug,'suppress_output'); fprintf('\n\n'); end;
            end
            
%             if ~strcmp(Debug, 'Bright') && ~strcmp(Debug, 'Quiet')
                obj = obj.norm_soundfield();
                obj = obj.extractSave_Zone('Bright');
                obj = obj.extractSave_Zone('Quiet');
                obj = obj.ErrorInBright_Calc();
                obj = obj.ErrorInQuiet_Calc();
                obj = obj.Contrast_Calc();
                obj.Angle_Bright_Field = Angle_Calc(obj, obj.Bright_Field);
                obj.Angle_Quiet_Field = Angle_Calc(obj, obj.Quiet_Field);
%             end
        end        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Soundfield Functions
        function obj = calc_Alpha_Coeffs(obj, Debug)
             if nargin < 2
                 Debug = '';
             end
             
            % Set up Orthogonal Basis Expansion for the Zones Here
            % Our goal is to find P so that we can apply it to the alpha coeffs                        
            obj = obj.w_build(Debug);  % Secondly we will set our weight function
            
            obj = obj.F_build(Debug);  % As a good starting point we then build our set of planewaves                        
            
            obj = obj.RG_build; % After which we perform orthonormalization with a weighted Gram-Schmidt process on the set of planewaves                     
                                % This gives us our upper triangular matrix and our orthonormal basis functions
                        
            obj = obj.C_build;  % Just before we find P we need C which is integral(Soundfield_Desired * Conjugate(G) * w)  
                                
            obj = obj.P_build(Debug);  % Finally we build our coefficient set P for the jth planewave function
            
            if contains(Debug,'DEBUG')
                [h,w] = size(obj.Soundfield_desired);
                F_r = reshape(obj.F, [h*w obj.N]);
                S_r = (F_r * obj.P');
                S = reshape( S_r, [h w 1]);
                obj = obj.createEmptySoundfield;
                obj.Soundfield_desired = S;
            end
            
            M = obj.getGlobalModeLimit();
            Alpha = zeros(1,2*M+1);
            
            for m = (-M:M)
                alpha = 1i^m * exp( -1i * m * obj.F_angles);

                Alpha(m+M+1) = sum( obj.P .* alpha );
            end
            obj.Alpha_Coeffs = Alpha;

        end

        function M0 = getGlobalModeLimit(obj)
            M0 = ceil(obj.k_global * obj.Radius * obj.ModeLimitFactor); 
        end
        
        function radius_p = getRadius_FromZones(obj)
           SpatialZones = [obj.Quiet_Zone obj.Bright_Zone];
           Q = length(SpatialZones); % Total number of zones
           
           % Find the radius of the smallest circle that encloses all spatial zones of interest.
           radius_p = SpatialZones(1).Radius_q + SpatialZones(1).getOriginInPolarCoords();
            for q = 2:Q 
                if (SpatialZones(q).Radius_q + SpatialZones(q).getOriginInPolarCoords()) > radius_p
                   radius_p = SpatialZones(q).Radius_q + SpatialZones(q).getOriginInPolarCoords();
                end
            end
        end                
     
%% Orthogonal Basis Expansion method functions
% Weighting function
        function obj = w_build(obj, Debug)
            if nargin < 2
                Debug = '';
            end
            [height, width] = size(obj.Soundfield_desired);
            obj.w  = zeros(height,width);
            obj.Sd = zeros(height,width);
            
%             if strcmp(Debug, 'Bright')
%                 obj.w = obj.w + obj.BrightZ_Weight;
%                 xx = (length(obj.Bright_Zone.Soundfield_d)/2-(width-1)/2):(length(obj.Bright_Zone.Soundfield_d)/2+(width-1)/2);
%                 yy = xx;
%                 obj.Sd = obj.Bright_Zone.Soundfield_d(yy,xx);
%             else
                obj.w = obj.w + obj.UnattendedZ_Weight; % Apply unattended zone weighting
                
                SpatialZones = [obj.Quiet_Zone obj.Bright_Zone];
                for q = 1:length(SpatialZones) % For each spatial zone
                    
                    %                     znSz = int16( SpatialZones(q).ZoneSize * obj.res );
                    %                     rpSz = int16( obj.ReproRegionSize * obj.res );
                    %                     x = ceil(SpatialZones(q).Origin_q.X * obj.res);
                    %                     y = ceil(SpatialZones(q).Origin_q.Y * obj.res);
                    %                     yy = ((-znSz(1)/2+1):(znSz(1)/2)) + y + rpSz(1)/2;
                    %                     xx = ((-znSz(2)/2+1):(znSz(2)/2)) + x + rpSz(2)/2;
                    
                    [xx,yy] = obj.getZoneIndices( SpatialZones(q) );
                    
                    if strcmp(SpatialZones(q).SourceType, 'quiet')
                        obj.w(yy,xx) =  SpatialZones(q).Soundfield_d_mask * (obj.QuietZ_Weight - obj.UnattendedZ_Weight);  % Apply quiet zone weighting
                    else
                        obj.w(yy,xx) = SpatialZones(q).Soundfield_d_mask * (obj.BrightZ_Weight - obj.UnattendedZ_Weight);   % Apply bright zone weighting
                    end
                    obj.w(yy,xx) = obj.w(yy,xx) + obj.UnattendedZ_Weight;
                    
                    
%                     znSz = ( SpatialZones(q).ZoneSize * obj.res );
%                     y = ((-znSz(1)/2+1):(znSz(1)/2));
%                     x = ((-znSz(2)/2+1):(znSz(2)/2));
%                     [X,Y] = meshgrid(x/obj.res,y/obj.res);
%                     W = abs( complex(X,Y) );
%                     W(W<min(SpatialZones(q).ZoneSize/2)) = min(SpatialZones(q).ZoneSize/2);
%                     W = W - min(W(:));
%                     W = W / max(W(:));
%                     obj.w = W.*50 + 1;
                    
                    obj.Sd(yy,xx) = SpatialZones(q).Soundfield_d .* SpatialZones(q).Soundfield_d_mask;

                end
%             end
        end
        
% Planewave Coefficient set 
        function obj = P_build(obj, Debug)
            if nargin < 2
                Debug = '';
            end
            obj.P = zeros(1,obj.N);
            if( rcond(obj.R) >= 1e-12 )
                obj.P = obj.C / obj.R';
            elseif contains(Debug,'IllCond')
                obj.P = obj.C / obj.R';
            end
        end
        
% Expansion Coefficient set
        function obj = C_build(obj)                  
            [height, width] = size(obj.Soundfield_desired);
            obj.C = zeros(1, obj.N);
                        
            for n = 1:obj.N
                G_field = reshape(obj.G(:,n), [height width]);
                field = obj.Sd .* G_field .* obj.w;
                obj.C(n) = sum(field(:));
            end
        end        

% Planewave set
        function obj = F_build(obj, Debug)
%             if strcmp(Debug, 'Bright')
%                 width=length(obj.Soundfield_desired);
%                 x = ((1:width) - width/2) / obj.res + obj.Bright_Zone.Origin_q.X;
%                 y = ((1:width) - width/2) / obj.res + obj.Bright_Zone.Origin_q.Y;
%                 [xx,yy] = meshgrid(x,y);
%                 X = complex(xx,yy);
%             else
                [height, width] = size(obj.Soundfield_desired);

%             end
            
            obj.F = zeros(height, width, obj.N);
            k = obj.k_global;
            
            %Create alternating basis planewaves centred around desired direction
            %             pos = 1:((obj.N - 1) / 2);
            %             neg = -1:-1:-((obj.N - 1) / 2);
            %             ind = [pos; neg];
            %             ind = [0 ind(:)' obj.N/2];
            %             obj.F_angles = ind * obj.Delta_phi + obj.Bright_Zone.SourceOrigin.Angle/180*pi;
            
            obj.F_angles = ((1:obj.N)-1) * obj.Delta_phi;
            
            x = ((1:width) - width/2 - 1) / obj.res;
            y = ((1:height) - height/2 - 1) / obj.res;
            [xx,yy] = meshgrid(x,y);
            X = complex(xx,yy);
            
            for n = 1:obj.N % For each planewave
                phi = obj.F_angles(n);
                
                Fn = exp( 1i * k * (cos(phi)*real(X) + sin(phi)*imag(X))) .* obj.Soundfield_desired_mask; %Planewave formula ( e^(i*(kx+ky+kz)) )
                
                Fn=Fn.*exp(-1i*angle(Fn(floor(end/2),floor(end/2)))); %norm to centre
                
                obj.F(:,:,n) = repmat(Fn, [1 1 1]);
            end
                
        end        

% Basis Functions with QR Factorisation & Orthogonalisation using Weighted Gram-Schmidt method
        function obj = RG_build(obj)
            [height, width] = size(obj.Soundfield_desired);
            obj.G = zeros(width * width, obj.N); 
            obj.R = zeros(obj.N, obj.N); 
                   
            F_all = reshape(obj.F, [width*height obj.N]);
            [obj.G, obj.R] = Orthogonal_Basis_Expansion.Gram_Schmidt(F_all, obj.w(:));            
           %[obj.G, obj.R] = Orthogonal_Basis_Expansion.Gram_Schmidt_mex(F_all, obj.w(:));
        end        
        
%% One-Dimensional sampling methods
        function obj = saveBrightZoneSample(obj)
            % First we want to extract the bright zone as a separate matrix
            
            
            % Second we want to rotate the matrix by the angle of the planewave
            
            
            % Thirdly we want to find the average planewave across the matrix


        end
        
        
%% Auxilary Methods
        function obj = addSpatialZone(obj, spatialZone, Distance, Angle_degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Check for overlapping soundfields here before adding it %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            SpatialZones = [obj.Quiet_Zone obj.Bright_Zone];
            
            if length(SpatialZones) < 1
                obj.k_global = spatialZone.getWavenumber();
                obj.res      = spatialZone.res;
            elseif obj.k_global ~= spatialZone.getWavenumber()
                error(['A spatial zone does not match the frequency/wavenumber of the global soundfield it is being added to.' ...
                      '\n Global Soundfield Frequency = %0.2fHz\tGlobal Soundfield Wavenumber = %0.2f' ...
                      '\n Spatial Zone Frequency = %0.2fHz\t\tSpatial Zone Wavenumber = %0.2f'], ...
                      obj.getFrequency(), obj.k_global, spatialZone.Frequency, spatialZone.getWavenumber());
                return;
            elseif obj.res ~= spatialZone.res
                error(['A spatial zone does not match the sampling resolution of the global soundfield it is being added to.' ...
                      '\n Global Soundfield Resolution = %0.2d samples/m' ...
                      '\n Spatial Zone Resolution = %0.2d samples/m'], ...
                      obj.res, spatialZone.res);
                return;
            end
            
            Angle_rads = Angle_degrees * pi / 180 ;
            [X, Y] = pol2cart(Angle_rads, Distance);
            spatialZone = spatialZone.setOrigin(X, Y);
            
            if strcmp(spatialZone.SourceType, 'quiet')
                obj.Quiet_Zone = spatialZone;
            elseif strcmp(spatialZone.SourceType, 'pw')
                obj.Bright_Zone = spatialZone;
            elseif strcmp(spatialZone.SourceType, 'ps')
                obj.Bright_Zone = spatialZone;
            elseif strcmp(spatialZone.SourceType, 'noise')
                obj.Bright_Zone = spatialZone;
            else
                error('Incorrect zone type');
            end
            
        end     
        
        function f = getFrequency(obj)
            f = obj.k_global * 343 / (2 * pi);
        end
        
        function obj = setReproRegionSize(obj, ReproRegionSize)
           obj.ReproRegionSize = ReproRegionSize;
        end
        
        function obj = setWavenumberFromChildZone(obj)
           obj.k_global = obj.Bright_Zone.getWavenumber();
        end
        
        function obj = setN(obj, N)
            if N>=1
                obj.N = N;
            else                
                obj = obj.setWavenumberFromChildZone(); %Incase the child spatial zones have been changed
                obj.N = 2 * obj.getGlobalModeLimit() + 1;
            end                
            obj.Delta_phi = 2 * pi / obj.N;
        end
        
        function [x,y] = getZoneIndices( obj, Zone )
            
            znSz = int16( Zone.ZoneSize * obj.res );
            rpSz = int16( size(obj.Soundfield_desired) );
            x_ = ceil( Zone.Origin_q.X * obj.res );
            y_ = ceil( Zone.Origin_q.Y * obj.res );
            y = ((-znSz(1)/2+1):(znSz(1)/2)) + y_ + rpSz(1)/2;
            x = ((-znSz(2)/2+1):(znSz(2)/2)) + x_ + rpSz(2)/2;
        
        end
        
        function obj = calcDesiredMask(obj)
            if strcmpi( obj.Geometry, 'circle' )
                width = int16(obj.res * obj.Radius * 2);
                xyvec = linspace(-obj.Radius, obj.Radius , width);
                [xx,yy] = meshgrid(xyvec,xyvec);
                obj.Soundfield_desired_mask = xx.^2 + yy.^2 <= (obj.Radius)^2 * ones(width, width);
            elseif contains( lower(obj.Geometry), 'rect' )
                width  = obj.ReproRegionSize(1) * obj.res;
                height = obj.ReproRegionSize(2) * obj.res;
                obj.Soundfield_desired_mask = ones(width, height);
            end
        end
        
        function obj = norm_soundfield(obj)
            zone_ = obj.Bright_Zone;
%             r = int16( zone_.Radius_q * obj.res );
%             r0= obj.Radius * obj.res;
%             x = ceil(zone_.Origin_q.X * obj.res);
%             y = ceil(zone_.Origin_q.Y * obj.res);
%             % xx = (-r:r-1) + x + r0;
%             % yy = (-r:r-1) + y + r0;
%             xx = (-r+1:r) + x + r0;
%             yy = (-r+1:r) + y + r0;
            [xx,yy] = obj.getZoneIndices( zone_ );
            f_ = obj.Soundfield_desired( yy, xx ) .* zone_.Soundfield_d_mask;
            
            obj.Soundfield_desired = obj.Soundfield_desired ./ max(abs(real(f_(:))));
        end
        
        function obj = extractSave_Zone(obj, zone)
            if strcmp(zone, 'Bright')
                zone_ = obj.Bright_Zone;
            elseif strcmp(zone, 'Quiet')
                zone_ = obj.Quiet_Zone;
            else
                return
            end
            
%             r = int16( zone_.Radius_q * obj.res );
%             r0= obj.Radius * obj.res;
%             x = ceil(zone_.Origin_q.X * obj.res);
%             y = ceil(zone_.Origin_q.Y * obj.res);
%             % xx = (-r:r-1) + x + r0;
%             % yy = (-r:r-1) + y + r0;
%             xx = (-r+1:r) + x + r0;
%             yy = (-r+1:r) + y + r0;
            
            [xx,yy] = obj.getZoneIndices( zone_ );
            
            result = obj.Soundfield_desired( yy, xx ) .* zone_.Soundfield_d_mask;

            if strcmp(zone, 'Bright')
                obj.Bright_Field = result;
            elseif strcmp(zone, 'Quiet')
                obj.Quiet_Field = result;
            end            
        end       
        
        function direction = Angle_Calc(obj, field)
            angles = angle( field );
            [U,V,W] = surfnorm( angles );
            W(W~=0) = 1;
            U = U .* W;
            V = V .* W;
            angles = complex( U, V );
            angles( abs(angles) > max( abs( angles(:) )/2 ) ) = 0;
            direction = angle( mean2( angles(angles~=0) ) ) / pi * 180;
        end        
        
        function obj = ErrorInBright_Calc(obj)

            desired = obj.Bright_Zone.Soundfield_d .* obj.Bright_Zone.Soundfield_d_mask;
            actual  = obj.Bright_Field;
            obj.MSE_Bright_Field = mean( abs(desired(:) - actual(:)) .^ 2 ) / mean( abs(desired(:)) .^ 2 );
            % We calculate the error for every degree of phase and average
            % the result
            %obj.MSE_Bright_Field = 0;
%             for p = 1:360
%                 desired = real((obj.Bright_Zone.Soundfield_d * exp(-1i*p/180*pi)) .* obj.Bright_Zone.Soundfield_d_mask);
%                 actual  = real(obj.Bright_Field * exp(1i*p/180*pi));
%                 
%                 % MSE Value
%                 obj.MSE_Bright_Field = obj.MSE_Bright_Field + ...
%                                        sum(sum( abs(desired - actual) .^ 2 )) / ...
%                                        sum(sum( abs(   desired      ) .^ 2 )) / 360;
%             end

            % Error in Decibels
            obj.Err_dB_Bright_Field = mag2db( obj.MSE_Bright_Field );
            
            % MSE Field
            obj.Bright_Error_Field =  (abs(desired - actual).^2) .* obj.Bright_Zone.Soundfield_d_mask./ ...
                                            abs(desired).^2 ;
        end
        
        function obj = ErrorInQuiet_Calc(obj)
            % We calculate the error for every degree of phase and average
            % the result
            quiet = zeros(size(obj.Quiet_Field));
            obj.Quiet_Error_Field = zeros(size(obj.Quiet_Field));
            obj.MSE_Quiet_Field = 0;
            
            for p = 1:360
                quiet = obj.Quiet_Field * exp(1i*p/180*pi);
                % MSE Field
                obj.Quiet_Error_Field = abs(real(quiet )).^2 + 1i * abs(imag(quiet )).^2;
                obj.Quiet_Error_Field = obj.Quiet_Error_Field ./ nnz(obj.Quiet_Zone.Soundfield_d_mask);
                
                % MSE Value
                obj.MSE_Quiet_Field = obj.MSE_Quiet_Field + ...
                                     (sum( real(obj.Quiet_Error_Field(:)) ) / 360 + ...
                                      sum( imag(obj.Quiet_Error_Field(:)) ) / 360) / 2;
                
            end
            
            % Error in Decibels
            obj.Err_dB_Quiet_Field = mag2db( obj.MSE_Quiet_Field );
            
            % MSE Field (remove outside information)
            obj.Quiet_Error_Field = obj.Quiet_Error_Field ./ obj.Quiet_Zone.Soundfield_d_mask;
        end
        
        function obj = Contrast_Calc(obj)
            % We calculate the contrast for every degree of phase and average
            % the result            
            bright = zeros(size(obj.Bright_Field));
            quiet  = zeros(size(obj.Quiet_Field));
            numer  = zeros(size(obj.Bright_Field));
            denom  = zeros(size(obj.Quiet_Field));
            obj.Acoustic_Brightness_Contrast = 0;
            
            for p = 1:360
                bright = obj.Bright_Field * exp(1i*p/180*pi);
                quiet  = obj.Quiet_Field  * exp(1i*p/180*pi);
                
                % MSE difference
                numer = abs(real(bright)).^2 + 1i * abs(imag(bright)).^2;
                denom = abs(real(quiet )).^2 + 1i * abs(imag(quiet )).^2;
                numer = (sum( real(numer(:)) ) + sum( imag(numer(:)) ))/2 / nnz(obj.Bright_Zone.Soundfield_d_mask);
                denom = (sum( real(denom(:)) ) + sum( imag(denom(:)) ))/2 / nnz(obj.Quiet_Zone.Soundfield_d_mask );
                
                % Acoustic Brightness Contrast Value
                obj.Acoustic_Brightness_Contrast = obj.Acoustic_Brightness_Contrast + numer / denom / 360;
            end
            
            % Acoustic Brightness Contrast in Decibels
            obj.Acoustic_Brightness_Contrast_dB = mag2db( obj.Acoustic_Brightness_Contrast );
                                    
        end
                
%% Plotting Functions
        function h = plotSoundfield(obj, field, colour, plotZones)
            if nargin < 4; plotZones = true; end;
            if nargin < 3; colour = parula; end;
            if nargin < 2; field = obj.Soundfield_desired; end;
            
            ZaxisSize = max(abs(real(field(:))));
            
            CaxisSize = ZaxisSize;
            if ~isempty(obj.Bright_Field(:))
                CaxisSize = max(abs(real(obj.Bright_Field(:))))*4/3;
            end
            
            %XYTick = [1 length(field)/4 length(field)/2 length(field)*3/4 length(field)];
            %XYTickLabel = { -length(field)/2/obj.res, -length(field)/4/obj.res, 0, length(field)/4/obj.res, length(field)/2/obj.res};
            XYTick = [1  length(field)/2  length(field)];
            XYTickLabel = { -length(field)/2/obj.res,  0,  length(field)/2/obj.res};
            
            %             subplot(1,2,1);
            h = surf(real(field),'EdgeColor','None');
            colormap(colour);
            view(2);
            title('Real',...
                'FontWeight','bold',...
                'FontSize',24,...
                'FontName','Arial');
            axis('square');
            box on;
            axis([1 length(field) 1 length(field)]);
            set(gca, 'XTick', XYTick); set(gca, 'XTickLabel', XYTickLabel);
            set(gca, 'YTick', XYTick); set(gca, 'YTickLabel', XYTickLabel);
            set(gca, 'FontSize', 16);
            set(gca, 'FontName', 'Arial');
            % Create xlabel
            xlabel('Width (metres)','FontSize',20,'FontName','Arial');
            % Create ylabel
            ylabel('Length (metres)','FontSize',20,'FontName','Arial');
            % Create zlabel
            zlabel('Amplitude','Visible','off');
            caxis([-CaxisSize CaxisSize]);
            if (strcmp(version('-release'), '2012b'))
                tick_str = 'Ztick';
            else                
                tick_str = 'Ticks';
            end
            colorbar(tick_str,[-1 -0.5 0 0.5 1],'FontSize',16, 'FontName','Arial');
           
            if plotZones
                obj.zonePlots(real(field));
            end
            %             subplot(1,2,2);
            %             surf(imag(field),'EdgeColor','None');
            %             view(2);
            %             title('Imag');
            %             axis('square');
            %             axis([1 length(field) 1 length(field)]);
            %             set(gca, 'XTick', XYTick); set(gca, 'XTickLabel', XYTickLabel);
            %             set(gca, 'YTick', XYTick); set(gca, 'YTickLabel', XYTickLabel);
            %             caxis([-CaxisSize CaxisSize]);
            %             colorbar;
            %             colormap(colour);
            %             obj.zonePlots();
            
        end
        
        function zonePlots(obj, field)
            if nargin < 2; field = obj.Soundfield_desired; end;
                
            % For a Global Zone outline plot
            GlobalZone_temp = obj.Bright_Zone;
            GlobalZone_temp.Origin_q = struct('X',0,'Y',0);
            GlobalZone_temp.Radius_q = obj.Radius;

            maxZ = max(real(field(:)));
            hold on
            th = 0:pi/obj.res:2*pi;
            SpatialZones = [obj.Quiet_Zone obj.Bright_Zone GlobalZone_temp];
            
            for q = 1:length(SpatialZones)
                Ox = SpatialZones(q).Origin_q.X * obj.res + length(field) / 2 + 1;
                Oy = SpatialZones(q).Origin_q.Y * obj.res + length(field) / 2 + 1;
                xunit = SpatialZones(q).Radius_q * obj.res * cos(th) + Ox;
                yunit = SpatialZones(q).Radius_q * obj.res * sin(th) + Oy;
                plot3(xunit, yunit, ones(length(xunit),1)*maxZ, 'Color', [0 0 1], 'LineWidth', 3);
            end            
            
            hold off
        end
        
    end
    
end

