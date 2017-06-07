classdef loudspeaker_setup
    % LOUDSPEAKER_SETUP - This class provides practical setup information
    %                     from a compatible desired multizone soundfield.
    %
    %
    %   Author: Jacob Donley, University of Wollongong, Australia
    %   Email: Jacob.Donley089@uowmail.edu.au
    %
    
    properties
        % Settings
        res = 50;                           % samples per metre %Resolution of soundfield
        Dimensionality = 2;                 % 2D, 2.5D, 3D
        Loudspeaker_Count = 32;             % Number of speakers used in the reproduction (Number of speakers required = ceil( Angular_Window * (2*M + 1) / (2*pi) ) ).
        Origin = [0,0];                     % [y,x]
        RoomSize = [];
        Radius = 1.5;                       % Radius of speaker layout
        Speaker_Arc_Angle = 180;            % Angle of the arc of speakers (360 = full circle, 180 = semi-circle)
        Angle_FirstSpeaker = 0;             % The angle where the first speaker of the array occurs
        Speaker_Array_Centre = 180;         % Angle of the centre of the loudspeaker arc
        Speaker_Array_Length = 0;           % Length of array (inclusive of end points)
        Speaker_Spacing = 0.01;             % The spacing between each consecutive loudspeaker
        Speaker_Array_Type = 'circle';      % 'circle' or 'line' or 'coprime'
        ExtendedField = false;
        Loudspeaker_Dimensions;
        Loudspeaker_Type = 'Genelec 8010A'; % 'Genelec 8010A' or 'Genelec 8020C' or 'Meyer MM-4XP' or 'Parametric'
        Loudspeaker_Object = [];
        %         DipoleDistance = 343/12/1000; % a distance of [343m/s / (12 x 1kHz)] will give approximately equal magnitude on one side of dipole array at 1kHz
        DipoleDistance = 343/(2*pi*2000);   % Follows formula in Donley et al, "Active Speech Control using Wave-Domain Processing with a Linear Wall of Dipole Secondary Sources", ICASSP, IEEE, 2017.
        k_global = 2000 /343*2*pi;          % Frequency in wavenumber
        c = 343;
        Multizone_Soundfield;               % multizone_soundfield object for reproduction
        %N_zone_samples = 1;                % The number of samples in each zone
        %Sampling_method = 'Centre';        % The layout of the samples in the zones. ('Centre', 'Random', 'Square', 'Circlular', 'All')
        
        % Results
        % 3D
        Loudspeaker_Weights = [];
        Loudspeaker_Locations = [];
        Loudspeaker_Directions = [];
        Soundfield_reproduced = [];         % The complex values of the reproduced sound field that is produced from this class.
        Soundfield_virtual = [];            % The complete virtual source soundfield
        Desired_Mask = [];
        %         Soundfield_Error = [];         % The error in the bright
        %         field
        % 2D
        %         Sound_sample_Avg = [];  %The average sound wave across the bright field in the direction of the desired planewave propagation
        %         Binaural_sample = [];   %TODO: Create binaural sample
        % 1D
        %         MSE_Bright_Field = 0;
        %         Err_dB_Bright_Field = 0;
        %         MSE_Quiet_Field = 0;
        %         Err_dB_Quiet_Field = 0;
        %         Acoustic_Brightness_Contrast = 0;
        %         Acoustic_Brightness_Contrast_dB = 0;
        Bright_Samples = [];
        Bright_Samples_Locations = [];
        Quiet_Samples  = [];
        Quiet_Samples_Locations = [];
        Bright_Sample = 0;
        Quiet_Sample  = 0;
        
        Acoustic_Contrast = 0;
        MSE_Bright = 0;
        Attenuation_dB = 0; %Log-normal mean and std error
    end
    
    properties (Access = private)
        
    end
    
    %% Contructor
    methods
        function obj = loudspeaker_setup(~)
            
        end
    end
    
    %% Private Methods
    methods (Access = private)
        
        function obj = createEmptySoundfield(obj)
            [width, height] = obj.getFieldSize();
            obj.Soundfield_reproduced = zeros(height,width);
        end
        
        function [width, height] = getFieldSize(obj)
            thisRadius = obj.Radius;
            if obj.Loudspeaker_Count == 1
                thisRadius = 0;
            end
            if isempty(obj.RoomSize)
                R = max([thisRadius, obj.Multizone_Soundfield.Radius]);
            end
            if obj.ExtendedField && isempty(obj.RoomSize)
                if ~isempty(strfind(obj.Speaker_Array_Type,'line')) % Mirror accross linear array instead
                    L = obj.Loudspeaker_Count;
                    if strcmpi(obj.Speaker_Array_Type, '2line')
                        L = L/2;
                    end
                    R = max([thisRadius, obj.Multizone_Soundfield.Radius]);                    
%                     R = max([R, (obj.Loudspeaker_Dimensions(1)*L)/2]);
                    width = int16(obj.res * 2*R * 2*1.25);
                    height = int16(obj.res * R * 2*1.5);
                else
                    width = int16(obj.res * 3*R * 2);
                    height = width;
                end
            elseif isempty(obj.RoomSize)
                width = int16(obj.res * R * 2);
                height = width;
            else
                height = int16(obj.RoomSize(1)*obj.res);
                width  = int16(obj.RoomSize(2)*obj.res);
            end
        end
        
    end
    
    %% Public Methods
    methods
        
        %%%%%%%%%%%%%%%%%%  Main Soundfield Creation Function  %%%%%%%%%%%%%%%%%%%%
        
        function obj = reproduceSoundfield(obj, Debug, Radius)
            %  if ~strcmp(Debug, 'SAMPLES_ONLY')
            if nargin < 2
                obj = obj.setRadius(obj.Radius);
                Debug = '';
            elseif nargin < 3 || isempty(Radius)
                obj = obj.setRadius(obj.Radius);
            else
                obj = obj.setRadius(Radius);
            end
            
            obj = obj.setWavenumberFromChild();
            
            obj = obj.createEmptySoundfield;
            
            obj = obj.save_Bright_Samples();
            obj = obj.save_Quiet_Samples();
            
            if ~strcmp(Debug, 'SAMPLES_ONLY')
                
                O = size(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
                x = ((1:O(2)*2) - O(2)) / obj.res - obj.Origin(2);
                y = ((1:O(1)*2) - O(1)) / obj.res - obj.Origin(1);
                [xx,yy] = meshgrid(x,y);
                
                obj.Soundfield_reproduced = obj.computeField(xx,yy);
                
            end
            
            % Equalise frequency dependent decay
            T = abs(besselh(0, obj.k_global)); % Magnitude compensation
            obj.Loudspeaker_Weights = obj.Loudspeaker_Weights .* T; %Remove effect of frequency dependent decay (2D ATF compensation?)
            
            % Normalise to the maximum of the absolute bright zone values
            bMask = obj.Multizone_Soundfield.Bright_Zone.Soundfield_d_mean_mask;
            BrightSamples = obj.Bright_Samples;
            BrightSamples(isnan(BrightSamples))=0;
            maxBrightVal = mean(abs(BrightSamples(bMask(:)))); % Changed to find the average as our goal is to have the bright zone fit on average.
            % Normalise phase about desired bright zone centre
            %             O = floor(size(obj.Bright_Samples)/2);
            %             centBrightAngle = angle(obj.Bright_Samples(O(1),O(2))) ...
            %                 -angle(obj.Multizone_Soundfield.Bright_Zone.Soundfield_d(O(1),O(2)));
            %             spkrAngleNorm = angle(obj.Loudspeaker_Weights(ceil(end/2)));%Normalise to centre of loudspeaker array
            %             adjAngle = 0;%spkrAngleNorm;% + centBrightAngle;
            
            % Normalise and align phase to bright zone
            obj.Loudspeaker_Weights     = obj.Loudspeaker_Weights   / maxBrightVal;% * exp(-1i*adjAngle); %Added to keep speaker weights tied with their response
            obj.Soundfield_reproduced   = obj.Soundfield_reproduced / maxBrightVal;% * exp(-1i*adjAngle);
            obj.Bright_Samples          = obj.Bright_Samples        / maxBrightVal;% * exp(-1i*adjAngle);
            obj.Quiet_Samples           = obj.Quiet_Samples         / maxBrightVal;% * exp(-1i*adjAngle);
            obj.Bright_Sample           = obj.Bright_Sample         / maxBrightVal;
            obj.Quiet_Sample            = obj.Quiet_Sample          / maxBrightVal;
            
            % Compute measures
            obj.Acoustic_Contrast = obj.getAcoustic_Contrast();
            obj.MSE_Bright = obj.getBrightError();
            obj.Attenuation_dB = obj.getMaxPossibleAttenuation();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [S,H] = computeField(obj, xx, yy, freqs)
            if nargin < 4
                freqs=[];
            end
            th2 = obj.Loudspeaker_Locations(:,1);
            r2 = obj.Loudspeaker_Locations(:,2) ;
            L = obj.Loudspeaker_Weights;
            di = obj.Loudspeaker_Directions/pi*180;
            xwid = size(xx,2);
            ywid = size(xx,1);
            Nspkrs = obj.Loudspeaker_Count;
            [th1,r1] = cart2pol(xx,yy);
            
            % Create matrices
            r1_  = repmat(r1 ,1,1,Nspkrs);
            th1_ = repmat(th1,1,1,Nspkrs);
            xx_  = repmat(xx ,1,1,Nspkrs);
            yy_  = repmat(yy ,1,1,Nspkrs);
            r2_  = permute( repmat(r2     ,1,ywid,xwid), [2 3 1] );
            th2_ = permute( repmat(th2    ,1,ywid,xwid), [2 3 1] );
            L_   = permute( repmat(conj(L),1,ywid,xwid), [2 3 1] );
            Di   = permute( repmat(di     ,1,ywid,xwid), [2 3 1] );
            
            th_ = th1_-th2_;
            r = (r1_.^2 + r2_.^2 - 2*r1_.*r2_.*cos(th_)).^(1/2);
            [Lxx,Lyy] = pol2cart(th2_,r2_);
            th1 = atan2( yy_-Lyy, xx_-Lxx );
            
            % Loudspeaker transfer functions
            if strcmp(obj.Loudspeaker_Type, 'Parametric')
                H = obj.Loudspeaker_Object.convolutionalModel( th1, r, Di, freqs);
            else
                if obj.Dimensionality == 2
                    H = 1i/4 * besselh(0, obj.k_global * r ); % 2D
                elseif obj.Dimensionality == 3
                    % H = 1i/4 * sqrt(pi./(2*obj.k_global * r)) .* besselh(0 + 0.5, obj.k_global * r ); % 3D
                    % H = exp( -1i*obj.k_global*r ) ./ (4*pi*r); % 3D
                    H = exp( 1i*obj.k_global*r ) ./ (4*pi*r); % 3D
                end
            end
            
            if (ndims(L_) == ndims(H)) && all(size(L_) == size(H))
                S = sum( L_.*H , 3 );
            else
                S=[];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Loudspeaker Functions
        function obj = calc_Loudspeaker_Weights(obj, Debug)
            if nargin < 2
                Debug = '';
            end
            
            if isempty( obj.Loudspeaker_Locations )
                error('Please set the loudspeaker locations first.') ;
            end
            
            obj = obj.setWavenumberFromChild();
            
            Q0 = obj.Loudspeaker_Count;
            
            if  Q0 > 1
                phi = obj.Speaker_Arc_Angle / 180 * pi;
                if strcmpi(obj.Speaker_Array_Type, '2line')
                    Q0_ = Q0/2;
                else
                    Q0_ = Q0;
                end
                delta_phi_s = phi / Q0_;
                M = obj.Multizone_Soundfield.getGlobalModeLimit;
                N = obj.Multizone_Soundfield.N;
                
                phi_p = permute(repmat(obj.Multizone_Soundfield.F_angles(:).',Q0,1,2*M+1),[1 3 2]);
                P_j = permute(repmat(obj.Multizone_Soundfield.P(:).',Q0,1,2*M+1),[1 3 2]);
                k=obj.k_global;
                m_ = repmat(-M:M,Q0,1,N);
                m=m_(:,:,1);
                phi_q = repmat(obj.Loudspeaker_Locations(:,1),1,2*M+1) + pi;
                R_q = repmat(obj.Loudspeaker_Locations(:,2),1,2*M+1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.Dimensionality == 2                    
                    % 2D
                    ATF = besselh(m, k * R_q ); % 2D
                elseif obj.Dimensionality == 3
                    % 3D
                    ATF = sqrt(pi./(2*k*R_q)) .* besselh(m + 0.5, k * R_q ); % 3D
                end
                
                % Expansion
                obj.Loudspeaker_Weights = sum(2 * exp(1i*m.*phi_q) * delta_phi_s .*   sum( P_j .* 1i.^m_ .* exp( -1i * m_ .* phi_p ),3) ...
                    ./ (1i*pi* ATF )  ,2).';
                
                %                 A=sum( P_j .* 1i.^m_ .* exp( -1i * m_ .* phi_p ),3);A=A(1,:);
                %                 m=(-M:M);
                %                 phi_p=obj.Multizone_Soundfield.F_angles(:);
                %                 P_j = obj.Multizone_Soundfield.P(:);
                %                 phi_q=obj.Loudspeaker_Locations(:,1)+pi;
                %                 R_q = obj.Loudspeaker_Locations(:,2);
                %
                %                 B = 2 * exp(1i*phi_q*m) * delta_phi_s;
                %                 C =  P_j.' * exp( -1i * phi_p * m ) .* 1i.^m ;
                %                 D = (1i*pi*besselh(repmat(m,Q0,1), k * repmat(R_q,1,2*M+1)));
                %                 Z = B .* kron(C,ones(Q0,1)) ./ D;
                %                 obj.Loudspeaker_Weights = sum(2 * exp(1i*m.*phi_q) * delta_phi_s .* C   ...
                %                                                 ./ (1i*pi*besselh(m, k * R_q))  ,2).';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else
                obj.Loudspeaker_Weights = 1;
            end
            
            obj.Loudspeaker_Weights = obj.Loudspeaker_Weights' ;
            
            if strcmpi(obj.Speaker_Array_Type, '2line') %this is a dipole setup (infinitesimal spacing)
                % Here we perform a delay and sum beamforming opperation to
                % create better dipole directivity for the '2line' case
                d = obj.DipoleDistance;
                f = obj.getFrequency;
                Q = obj.Loudspeaker_Count;
                LW= obj.Loudspeaker_Weights;
                dpF = exp(-1i*pi/2) * ...
                    [exp( 1i*pi), ...
                    exp( 1i*pi*2*d*f/obj.c)];
                dpFs = reshape(repmat(dpF,Q/2,1),[],1);
                LW(1:end/2) = LW(end:-1:end/2+1);
                obj.Loudspeaker_Weights = LW .* dpFs;
            end
            
        end
        
        
        function obj = calc_Loudspeaker_Locations(obj)
            Q0 = obj.Loudspeaker_Count;
            Rl = obj.Radius;
            
            if strcmpi(obj.Speaker_Array_Type, 'circle')
                phi = obj.Speaker_Arc_Angle / 180 * pi;
                delta_phi_s = phi / Q0;
                first_spkr_offset = obj.Angle_FirstSpeaker / 180 * pi;
                phi_q = ((1:Q0) - 1) * delta_phi_s + pi + first_spkr_offset;
                Rl_ = repmat(Rl, length(phi_q), 1);
                obj.Loudspeaker_Locations = [(phi_q(:) - pi) Rl_];
                
            elseif strcmpi(obj.Speaker_Array_Type, 'line')
                wid = obj.Loudspeaker_Dimensions(1);
                space = obj.Speaker_Spacing;
                centre = obj.Speaker_Array_Centre;
                L = obj.Loudspeaker_Count;
                len = (wid+space)*L - space;
                x = repmat(obj.Radius, 1, L);
                
                if L == 1
                    a2first = obj.Angle_FirstSpeaker/180*pi;
                    a2center = centre/180*pi;
                    spkr_rho = x(1) / cos(a2first - a2center);
                    obj.Loudspeaker_Locations = [a2first spkr_rho];
                else
                    
                    if ~mod(L,2)
                        %                         y = ((wid+space)/2):(wid+space):len/2;
                        y = (1/2:L/2)*(wid+space);
                        y = [flip(-y), y];
                    else
                        y = (1:(L-1)/2)*(wid+space);
                        y = [flip(-y), 0, y];
                    end
                    
                    [spkr_theta, spkr_rho]=cart2pol(x,y);
                    spkr_theta = spkr_theta + centre/180*pi;
                    
                    obj.Loudspeaker_Locations = [spkr_theta(:), spkr_rho(:)];
                    
                    obj.Angle_FirstSpeaker = obj.Loudspeaker_Locations(1,1) / pi * 180;
                    obj.Speaker_Arc_Angle = obj.Loudspeaker_Locations(end,1) / pi * 180 - obj.Angle_FirstSpeaker;
                    
                    LLends = obj.Loudspeaker_Locations([1 end],:);
                    [LLx,LLy] = pol2cart(LLends(:,1),LLends(:,2));
                    obj.Speaker_Array_Length = sum(diff([LLx, LLy]).^2).^0.5; % inclusive of end points
                end
                
            elseif strcmpi(obj.Speaker_Array_Type, '2line') %this is a dipole setup (infinitesimal spacing)
                
                wid = obj.Loudspeaker_Dimensions(1);
                dist_=obj.DipoleDistance;
                space = obj.Speaker_Spacing;
                centre = obj.Speaker_Array_Centre;
                Ltot = obj.Loudspeaker_Count;
                if mod(Ltot,2)
                    error('Not an even number of loudspeakers. ''2line'' array type requires an even number.');
                end
                L = Ltot/2;
                len = (wid+space)*L - space;
                x = [repmat(obj.Radius+dist_/2, 1, L), ...
                    repmat(obj.Radius-dist_/2, 1, L)];
                
                if L == 1
                    a2first = obj.Angle_FirstSpeaker/180*pi;
                    a2center = centre/180*pi;
                    spkr_rho = x(1) / cos(a2first - a2center);
                    obj.Loudspeaker_Locations = [a2first spkr_rho];
                else
                    
                    y = ((wid+space)/2):(wid+space):len/2;
                    y = [flip(-y), y];
                    
                    y = [y, flip(y)];
                    
                    [spkr_theta, spkr_rho]=cart2pol(x,y);
                    spkr_theta = spkr_theta + centre/180*pi;
                    
                    obj.Loudspeaker_Locations = [spkr_theta(:), spkr_rho(:)];
                    
                    [AngFir,~] = cart2pol(mean(x([1 end])),y(1));
                    [AngLas,~] = cart2pol(mean(x(end/2 + [1 -1])),y(end/2));
                    
                    obj.Angle_FirstSpeaker = AngFir/pi*180 + centre;
                    obj.Speaker_Arc_Angle = AngLas/pi*180 + centre - obj.Angle_FirstSpeaker;
                    
                    obj.Speaker_Array_Length = sum(abs(y([1 end/2]))); % inclusive of end points
                end
                
            elseif strcmpi(obj.Speaker_Array_Type, 'coprime')
                L = obj.Loudspeaker_Count;
                
                if ~isprime(L+1) && ~isprime(L)
                    error('Number of loudspeakers is not one less than a prime number.');
                end
                
                wid = obj.Loudspeaker_Dimensions(1);
                space = obj.Speaker_Spacing;
                centre = obj.Speaker_Array_Centre;
                len = (wid+space)*L - space;
                x = repmat(obj.Radius, L, 1);
                
                ends_ = (len/2-wid/2);
                y1 =  linspace(ends_,-ends_,L/2+1);
                y1(end) = [];
                y2 =  linspace(ends_,-ends_,L/2+2);
                y_ = -y2(1);
                y2(1)=[];
                y2(end)=[];
                
                y = [y1; y2];
                if iseven(L)
                    y = y(:);
                else
                    y = [y(:);y_];
                end
                
                [spkr_theta, spkr_rho]=cart2pol(x,y);
                spkr_theta = spkr_theta + centre/180*pi;
                
                obj.Loudspeaker_Locations = [spkr_theta(:), spkr_rho(:)];
                
                obj.Angle_FirstSpeaker = obj.Loudspeaker_Locations(1,1) / pi * 180;
                obj.Speaker_Arc_Angle = obj.Loudspeaker_Locations(end,1) / pi * 180 - obj.Angle_FirstSpeaker;
                
            end
            
            % Determine the directions of the loudspeakers depending on the
            % type of array and the position of the speakers in that array.
            obj = obj.calc_Loudspeaker_Directions();
            
        end
        
        
        function obj = calc_Loudspeaker_Directions(obj, directions)
            if nargin >= 2
                obj.Loudspeaker_Directions = directions;
            else
                Q0 = obj.Loudspeaker_Count;
                
                if strcmp(obj.Speaker_Array_Type, 'circle')
                    obj.Loudspeaker_Directions = obj.Loudspeaker_Locations(:,1) - pi;
                    
                elseif strcmp(obj.Speaker_Array_Type, 'line') ...
                        || strcmp(obj.Speaker_Array_Type, 'coprime')
                    obj.Loudspeaker_Directions = ones(Q0,1) * (obj.Speaker_Array_Centre - 180) / 180 * pi;
                    
                end
            end
        end
        
        %% Auxilary Methods
        function obj = addMultizone_Soundfield(obj, multizone_soundfield )
            
            obj.Multizone_Soundfield = [];
            
            if length(obj.Multizone_Soundfield) < 1
                obj.res      = multizone_soundfield.res;
            elseif obj.res ~= multizone_soundfield.res
                error(['A multizone soundfield does not match the sampling resolution of the previously added soundfield(s).' ...
                    '\n Previous Soundfield Resolution = %0.2d samples/m' ...
                    '\n New Soundfield Resolution = %0.2d samples/m'], ...
                    obj.res, multizone_soundfield.res);
            end
            
            obj.Multizone_Soundfield = multizone_soundfield;
            
            obj.k_global = obj.Multizone_Soundfield.k_global;
            obj = obj.setRadius(obj.Multizone_Soundfield.Radius);
            
        end
        
        function obj = addLoudspeaker_Object(obj, loudspeaker_object)
            obj.Loudspeaker_Object  = loudspeaker_object;
            obj.Loudspeaker_Object.res = obj.res;
            obj = obj.createEmptySoundfield;
            obj.Loudspeaker_Object.field_size = size(obj.Soundfield_reproduced) / obj.res;
            if isfield(obj.Loudspeaker_Object,'set_fd')
                obj.Loudspeaker_Object = obj.Loudspeaker_Object.set_fd(obj.k_global*obj.c/2/pi);
            end
        end
        
        function obj = calcDesiredMask(obj)
            %             [w,h] = obj.getFieldSize();
            %             x = linspace(-obj.Radius, obj.Radius , w);
            %             y = linspace(-obj.Radius, obj.Radius , h);
            %             [xx,yy] = meshgrid(xyvec,xyvec);
            [O(2), O(1)] = obj.getFieldSize(); %Set the centre of the zone for indexing
            O = double(O)/2;
            x = ((1:O(2)*2) - O(2)) / obj.res - obj.Origin(2);
            y = ((1:O(1)*2) - O(1)) / obj.res - obj.Origin(1);
            [xx,yy] = meshgrid(x,y);
            obj.Desired_Mask = xx.^2 + yy.^2 <= (obj.Multizone_Soundfield.Radius)^2 * ones(numel(y), numel(x));
        end
        
        function f = getFrequency(obj)
            f = obj.k_global * obj.c / (2 * pi);
        end
        
        function obj = setWavenumberFromChild(obj)
            obj.k_global = obj.Multizone_Soundfield.k_global;
            if ~isempty(obj.Loudspeaker_Object) && isfield(obj.Loudspeaker_Object,'set_fd')
                obj.Loudspeaker_Object = obj.Loudspeaker_Object.set_fd(obj.k_global*obj.c/2/pi);
            end
        end
        
        function obj = setRoomSize(obj, RoomSize)
            obj.RoomSize = RoomSize;
        end
        
        function obj = setRadius(obj, Radius)
            obj.Radius = Radius;
            obj = obj.calcDesiredMask();
        end
        
        
        function obj = setLoudspeakerType(obj, Type)
            obj.Loudspeaker_Type = Type;
            
            % Dimensions: [width, depth, height]
            if strcmp(Type, 'Genelec 8010A')
                obj.Loudspeaker_Dimensions = [0.121, 0.116, 0.195]; %Genelec 8010A
                
            elseif strcmp(Type, 'Genelec 8020C')
                obj.Loudspeaker_Dimensions = [0.151, 0.142, 0.242]; %Genelec 8020C
                
            elseif strcmp(Type, 'Meyer MM-4XP')
                obj.Loudspeaker_Dimensions = [0.103, 0.145, 0.103]; %Meyer MM-4XP
                
            elseif strcmp(Type, 'Parametric')
                obj.Loudspeaker_Dimensions = obj.Loudspeaker_Object.Spkr_Dimensions; %Parametric Array Loudspeaker (PAL) Dimensions
                
            else
                error('Loudspeaker type not supported.');
            end
        end
        
        function obj = setLoudspeakerSpacing(obj, Speaker_Arr_Centre, LoudspeakerSpacing)
            width = obj.Loudspeaker_Dimensions(1);
            if nargin < 2
                centre = obj.Speaker_Array_Centre;
            else
                centre = Speaker_Arr_Centre;
            end
            if nargin < 3
                space = obj.Speaker_Spacing;
            else
                space = LoudspeakerSpacing;
            end
            
            circumf = obj.Loudspeaker_Count * (space + width);
            
            new_Speaker_Arc_Angle = circumf / obj.Radius * 180/pi;
            if new_Speaker_Arc_Angle > obj.Speaker_Arc_Angle && ~isnan(obj.Speaker_Arc_Angle)
                wrnCol = [255,100,0]/255;
                wrnTxt = ['The previously specified loudspeaker arc is too small to accommodate the desired number of loudspeakers.' 10, ...
                    'Either:' 10 ...
                    9	'1) Reduce the number of loudspeakers.' 10 ...
                    9  '2) Use smaller loudspeakers.  or' 10 ...
                    9  '3) Increase the size of the loudspeaker arc.' 10, ...
                    'The loudspeaker arc will be adjusted now...' 10, ...
                    9  'Old loudspeaker arc: ~' num2str(round(obj.Speaker_Arc_Angle,1)) '°' 10, ...
                    9  'New loudspeaker arc: ~' num2str(round(new_Speaker_Arc_Angle,1)) '°' 10 10];
                cprintf(wrnCol, wrnTxt);
            end
            obj.Speaker_Arc_Angle = new_Speaker_Arc_Angle;
            if obj.Loudspeaker_Count > 1
                obj.Angle_FirstSpeaker = centre - obj.Speaker_Arc_Angle/2 + (((space + width)/2) / obj.Radius * 180/pi );
            end
        end
        
        
        function obj = save_Bright_Samples(obj)
            [obj.Bright_Samples, ...
                obj.Bright_Samples_Locations, ...
                obj.Bright_Sample] ...
                = getZoneSamples( obj, obj.Multizone_Soundfield.Bright_Zone );
        end
        
        function obj = save_Quiet_Samples(obj)
            [obj.Quiet_Samples, ...
                obj.Quiet_Samples_Locations, ...
                obj.Quiet_Sample] ...
                = getZoneSamples( obj, obj.Multizone_Soundfield.Quiet_Zone );
        end
        
        function [samples,locations,meansample] = getZoneSamples(obj, zone)
            % % If number is not a square number change method to random
            % if ~~rem(obj.N_zone_samples, sqrt(obj.N_zone_samples))
            % obj.Sampling_method = 'Random';
            % else
            %
            % end
            %
            % if strcmp(obj.Sampling_method, 'Random')
            %
            % end
            
            O = zone.Radius_q * obj.res;
            znSz = zone.ZoneSize * obj.res;
            %Oz = size(obj.Soundfield_reproduced) / 2 / obj.res; %Set the centre of the zone for indexing
            
            y = (-znSz(1)/2+1:znSz(1)/2) / obj.res;
            x = (-znSz(2)/2+1:znSz(2)/2) / obj.res;
            mask = zone.Soundfield_d_mask;
            maskNaN = 1*mask;
            maskNaN(~mask)=nan;
            
            [xx,yy] = meshgrid(x,y);
            xx = (xx + zone.Origin_q.X ) .* maskNaN ;
            yy = (yy + zone.Origin_q.Y ) .* maskNaN;
            
            [locations(:,:,1), ...
                locations(:,:,2)] = cart2pol(xx,yy);
            samples = obj.computeField(xx,yy);
            meansample = mean( abs( samples(mask) ) );
            
        end
        
        function contrast = getAcoustic_Contrast(obj)
            bMask = obj.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
            qMask = obj.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask;
            contrast = mean( abs( obj.Bright_Samples(bMask) ) .^ 2 ) / mean( abs( obj.Quiet_Samples(qMask) ) .^ 2 );
        end
        
        function MSE = getBrightError(obj)
            
            [desired_vec, actual_vec] = obj.getNormalisedFieldVectors;
            MSE = sum( abs(desired_vec - actual_vec) .^ 2 ) / sum( abs(desired_vec) .^ 2 );
            
        end
        
        function Atten_dB = getMaxPossibleAttenuation(obj)
            [desired_vec, actual_vec] = obj.getNormalisedFieldVectors;
            
            att = abs( desired_vec - actual_vec );
            dbConv = 20/log(10);
            m = mean(att.');
            v = var(att.');
            mu_ = log( m ./ sqrt( 1 + v./m.^2 ) ) *dbConv;
            sigm = sqrt( log( 1 + v./m.^2 ) ) *dbConv;
            
            Atten_dB = [mu_, sigm]; %lognormal mean and std error
        end
        
        function [Sd, Sa] = getNormalisedFieldVectors(obj)
            zoneMask = obj.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
            O = floor(size(zoneMask)/2);
            
            desired = obj.Multizone_Soundfield.Bright_Zone.Soundfield_d .* zoneMask;
            actual  = obj.Bright_Samples .* zoneMask;
            
            desired_vec=desired(zoneMask~=0);
            actual_vec=actual(zoneMask~=0);
            
            desired_adj = 1 / rms(abs(desired_vec)) .* exp(-1i* angle(desired(O(1),O(2))));
            actual_adj  = 1 / rms(abs(actual_vec)) .* exp(-1i* angle(actual(O(1),O(2))));
            
            Sd = desired_vec .* desired_adj;
            Sa = actual_vec .* actual_adj;
        end
        
        %% Plotting Functions
        function h = plotSoundfield(obj, field, colour, realistic, details)
            if nargin < 5; details.DrawDetails = false; end
            if nargin < 4; realistic = false; end
            if nargin < 3; colour = 'default'; end
            if nargin < 2; field = obj.Soundfield_reproduced; end
            
            ZaxisSize = max(abs(real(field(~isinf(field(:))))));
            f = real(field) .* obj.Desired_Mask;
            if min(real(f(:))) < 0
                CaxisSize = [-1 1].*ZaxisSize;
                if abs(min(real(f(:)))+pi) < 1e-3
                    CaxisSize = CaxisSize./obj.res.*pi;
                end
            else
                field = mag2db(field);
                fb = mag2db(abs(obj.Bright_Samples));
                fq = mag2db(abs(obj.Quiet_Samples));
                CaxisSize(1) = floor(min(fq(~isnan(fq))));
                CaxisSize(2) = ceil(max(fb(~isnan(fb))));
            end
            
            if ~details.DrawDetails
                details.NTicks = [5 5];
            end
            
            lenDim = size(field,2);
            XTick = fix(linspace(-1,1,details.NTicks(1))*lenDim/obj.res)*obj.res/2 + lenDim/2;
            XTickLabel = (XTick-lenDim/2)/obj.res-obj.Origin(2);
            XTickLabel = num2cell(XTickLabel);
            if diff(details.NTicks)<0
                for i=1:2:numel(XTickLabel)
                    XTickLabel{i}=[];
                end
            end
            
            lenDim = size(field,1);
            YTick = fix(linspace(-1,1,details.NTicks(2))*lenDim/obj.res)/2*obj.res + lenDim/2;
            YTickLabel = (YTick-lenDim/2)/obj.res-obj.Origin(1);
            YTickLabel = num2cell(YTickLabel);
            
            h = surf(real(field),'EdgeColor','None');
            if realistic
                drawnow; pause(0.05);
                h.FaceAlpha = 0.5;
            end
            if contains(colour,'scientific')
                colormap(cmap( strrep(colour,'scientific_','') ));
            else
                colormap(colour);
            end
            view(2);
            axis equal;
            box on;
            if realistic
                room = details.SYS.Room_Setup.Room_Size([2 1])*obj.res;
                axis([(size(field,2)/2-room(1)/2) (size(field,2)/2+room(1)/2) (size(field,1)/2-room(2)/2) (size(field,1)/2+room(2)/2)]);
            else
                axis([0 size(field,2)+1 0 size(field,1)+1]);
            end
            set(gca, 'XTick', XTick); set(gca, 'XTickLabel', XTickLabel);
            set(gca, 'YTick', YTick); set(gca, 'YTickLabel', YTickLabel);
            set(gca, 'TickDir', 'out' );
            set(gca, 'FontSize', 14);
            set(gca, 'FontName', 'cambria');
            %Create Title
            title('Pressure Soundfield',...
                'FontWeight','bold',...
                'FontSize',16,...
                'FontName','Times');
            % Create xlabel
            xlabel('Width (m)','FontSize',16,'FontName','Times');
            % Create ylabel
            ylabel('Length (m)','FontSize',16,'FontName','Times');
            % Create zlabel
            zlabel('Amplitude','Visible','off');
            caxis(CaxisSize);
            
            hCB = colorbar;
            hCB.Ticks = interp1(1:length(caxis),caxis,linspace(1,length(caxis),5));
            if sum(CaxisSize) == 0 && abs(min(real(f(:)))+pi) >= 1e-3
                hCB.TickLabels = num2str(linspace(CaxisSize(1),CaxisSize(2),5)' ./ obj.res);
            else
                hCB.TickLabels = num2str(linspace(CaxisSize(1),CaxisSize(2),5)' );
            end
            hCB.FontName = 'Times';
            hCB.FontSize = 16;
            
            details.ZoneGeometry = obj.Multizone_Soundfield.Geometry;
            details.ReproRegionSize = obj.Multizone_Soundfield.ReproRegionSize;
            
            obj.zonePlots(real(field),realistic, details);
            obj.sourcePlots(real(field),realistic);
            if details.DrawDetails
                obj.detailPlots(real(field),realistic,details);
            end
            
        end
        
        function zonePlots(obj, field, realistic, details)
            if nargin < 4; details.DrawDetails = false; end
            if nargin < 3; realistic = false; end
            if nargin < 2; field = obj.Soundfield_reproduced; end
            
            if details.DrawDetails
                LineWid = details.zoneLineWid;
            else
                LineWid = 1;
            end
            
            % For a Global Zone outline plot
            GlobalZone_temp = obj.Multizone_Soundfield.Bright_Zone;
            GlobalZone_temp.Origin_q = struct('X',0,'Y',0);
            GlobalZone_temp.Radius_q = obj.Multizone_Soundfield.Radius;
            GlobalZone_temp.ZoneGeometry = obj.Multizone_Soundfield.Geometry;
            GlobalZone_temp.ZoneSize = obj.Multizone_Soundfield.ReproRegionSize;
            
            if realistic
                maxZ = 0;
            else
                maxZ = max(real(field(:)));
            end
            hold on
            th = 0:pi/obj.res:2*pi;
            SpatialZones = [obj.Multizone_Soundfield.Quiet_Zone ...
                obj.Multizone_Soundfield.Bright_Zone ...
                GlobalZone_temp ];
            linestyles = {'-';'-';'--'}; % {Quiet Zone, Bright Zone, Reproduction Region}
            for q = 1:length(SpatialZones)
                O = [flip(size(field))/2 0] + [flip(obj.Origin)*obj.res 0];
                CentreP=[SpatialZones(q).Origin_q.X, ...
                    SpatialZones(q).Origin_q.Y, maxZ];
                if contains( lower(SpatialZones(q).ZoneGeometry), 'circle' )
                    plotCircle(obj, O, ...
                        CentreP, ...
                        SpatialZones(q).Radius_q, ...
                        {linestyles{q}, 'Color', [0 0 0], 'LineWidth', LineWid});
                elseif contains( lower(SpatialZones(q).ZoneGeometry), 'rect' )
                    plotRectangle(obj, O, ...
                        CentreP, ...
                        SpatialZones(q).ZoneSize, ...
                        0, ... Direction (angle of rotation)
                        {linestyles{q}, 'Color', [0 0 0], 'LineWidth', LineWid});
                end
                if details.DrawDetails
                    xunit = LineWid*3*cos(th) + CentreP(1)*obj.res+O(1);
                    yunit = LineWid*3*sin(th) + CentreP(2)*obj.res+O(2);
                    patch(xunit, yunit, ones(length(xunit),1)*maxZ, [0 0 0], 'LineStyle', 'none');
                end
            end
            
            hold off
        end
        
        function plotCircle(obj, Origin, Centre, Radius, plotArgs)
            th = 0:pi/obj.res:2*pi; %full circle
            C(1) = Centre(1) * obj.res + Origin(1);
            C(2) = Centre(2) * obj.res + Origin(2);
            C(3) = Centre(3)           + Origin(3);
            R = Radius * obj.res;
            Xunits = R * cos(th) + C(1);
            Yunits = R * sin(th) + C(2);
            Zunits = Xunits*0 + C(3); %planar circle plot
            hold on;
            plot3(Xunits, Yunits, Zunits, plotArgs{:});
            hold off;
        end
        
        function plotRectangle(obj, Origin, Centre, Size, Angle, plotArgs)
            C(1) = Centre(1) * obj.res + Origin(1);
            C(2) = Centre(2) * obj.res + Origin(2);
            C(3) = Centre(3)           + Origin(3);
            px_=[-1,1,1,-1;-1,1,1,-1];
            py_=[-1,-1,1,1;-1,-1,1,1];
            pz_=[-1,-1,-1,-1;1,1,1,1];
            Xunits = px_*Size(1)/2*obj.res + C(1);
            Yunits = py_*Size(2)/2*obj.res + C(2);
            Zunits = Xunits*0 + C(3); %planar rectangle plot
            hold on;
            plot3(Xunits([1:end 1]), Yunits([1:end 1]), Zunits([1:end 1]), plotArgs{:});
            hold off;
        end
        
        function detailPlots(obj, field, realistic, details)
            if nargin < 4; details.DrawDetails = false; end
            if nargin < 3; realistic = false; end
            if nargin < 2; field = obj.Soundfield_reproduced; end
            if realistic
                maxZ = 0;
            else
                maxZ = max(real(field(:)));
            end
            hold on
            
            LineWid = details.arrowLineWid;
            arrowLen = details.arrowLength;
            arrowTipAng = details.arrowAngle;
            endBuf = details.arrowBuffer / obj.res;
            fontSZ = details.lblFontSize;
            
            %Origin
            O = [length(field) / 2; ...
                length(field) / 2; ...
                maxZ * obj.res];
            
            % Zone Centres
            [Bth,Br]=cart2pol(obj.Multizone_Soundfield.Bright_Zone.Origin_q.X, ...
                obj.Multizone_Soundfield.Bright_Zone.Origin_q.Y);
            [Bx,By] = pol2cart(Bth,Br);
            [Bx_,By_] = pol2cart(Bth,Br-endBuf*2);
            BZ_Centre = [ Bx, By, 0]' * obj.res + O;
            BZ_Centre_ = [ Bx_, By_, 0]' * obj.res + O;
            
            [Qth,Qr]=cart2pol(obj.Multizone_Soundfield.Quiet_Zone.Origin_q.X, ...
                obj.Multizone_Soundfield.Quiet_Zone.Origin_q.Y);
            [Qx,Qy] = pol2cart(Qth,Qr);
            [Qx_,Qy_] = pol2cart(Qth,Qr-endBuf*2);
            QZ_Centre = [ Qx, Qy, 0]' * obj.res + O;
            QZ_Centre_ = [ Qx_, Qy_, 0]' * obj.res + O;
            
            % Zone Edges
            BZ_Edge = [obj.Multizone_Soundfield.Bright_Zone.Origin_q.X; ...
                obj.Multizone_Soundfield.Bright_Zone.Origin_q.Y + ...
                (obj.Multizone_Soundfield.Bright_Zone.Radius_q-endBuf); ...
                0]* obj.res + O;
            QZ_Edge = [obj.Multizone_Soundfield.Quiet_Zone.Origin_q.X; ...
                obj.Multizone_Soundfield.Quiet_Zone.Origin_q.Y - ...
                (obj.Multizone_Soundfield.Quiet_Zone.Radius_q-endBuf); ...
                0]* obj.res + O;
            
            %Reproduction Edge
            Repro_Edge_Angle = -45;%degrees
            [Repro_Edge(1), ...
                Repro_Edge(2), ...
                Repro_Edge(3)] = sph2cart(Repro_Edge_Angle/180*pi,0,(obj.Multizone_Soundfield.Radius-endBuf)*obj.res);
            Repro_Edge = Repro_Edge' + O;
            
            %Loudspeaker Arc
            [SpkrArc_Edge(1), ...
                SpkrArc_Edge(2), ...
                SpkrArc_Edge(3)] = sph2cart(obj.Speaker_Array_Centre/180*pi,0,(obj.Radius-endBuf)*obj.res);
            SpkrArc_Edge = SpkrArc_Edge' + O;
            
            % Plot detailed coordinates
            arrowSpecs={'Length',arrowLen,'Width',LineWid,'TipAngle',arrowTipAng,'BaseAngle',(90+arrowTipAng)/2};
            arrow(O,BZ_Centre_,arrowSpecs{:});
            arrow(O,QZ_Centre_,arrowSpecs{:});
            arrow(BZ_Centre,BZ_Edge,arrowSpecs{:});
            arrow(QZ_Centre,QZ_Edge,arrowSpecs{:});
            arrow(O,Repro_Edge,arrowSpecs{:});
            arrow(O,SpkrArc_Edge,arrowSpecs{:});
            pwAngle=obj.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle;
            drawPlaneWaveArrow(obj,BZ_Centre,BZ_Edge,pwAngle,arrowSpecs); %planewave arrow
            
            % Label arrows
            addText2Arrow(obj,'\it{r_z}',O,BZ_Centre_,fontSZ, -10,0);
            addText2Arrow(obj,'\it{r_z}',O,QZ_Centre_,fontSZ,-10,0);
            addText2Arrow(obj,'\it{r}',BZ_Centre,BZ_Edge,fontSZ,10,0);
            addText2Arrow(obj,'\it{r}',QZ_Centre,QZ_Edge,fontSZ,-10,0);
            addText2Arrow(obj,'\it{R}',O,Repro_Edge,fontSZ,10,0);
            addText2Arrow(obj,'\it{R_l}',O,SpkrArc_Edge,fontSZ,-10,0);
            
            % Plot angles
            R = obj.Multizone_Soundfield.Radius * obj.res *0.15;
            drawAngle(obj,O,Bth,obj.Speaker_Array_Centre/180*pi,R,LineWid,maxZ)
            addText2Angle(obj,'\it{\phi}',O,Bth,obj.Speaker_Array_Centre/180*pi,fontSZ,R+5,maxZ)
            drawAngle(obj,O,Bth,Qth,R,LineWid,maxZ)
            addText2Angle(obj,'\it{\pi}',O,Bth,Qth,fontSZ,R+5,maxZ)
            drawAngle(obj,BZ_Centre,pwAngle,-Bth,R,LineWid,maxZ)
            addText2Angle(obj,'\it{\theta}',BZ_Centre,pwAngle,-Bth,fontSZ,R+5,maxZ)
            
            % Auxilary components
            %Alias Point
            Q = QZ_Centre;
            B = BZ_Centre;
            Rl = (obj.Radius);
            r1 = obj.Multizone_Soundfield.Bright_Zone.Radius_q;
            r2 = obj.Multizone_Soundfield.Quiet_Zone.Radius_q;
            L = obj.Loudspeaker_Count;
            phi_L=obj.Speaker_Arc_Angle/180*pi;
            Rprime = max([r1+Br,r2+Qr]);
            theta = pwAngle/180*pi;
            arrLen = obj.Speaker_Array_Length;
            k = obj.k_global;
            Bv = B-O;
            Qv = Q-O;
            
            rotMat_90 = [cos(pi/2) -sin(pi/2) 0;...
                sin(pi/2)  cos(pi/2) 0;...
                0           0      1];
            
            
            
            if strcmpi('line',obj.Speaker_Array_Type)
                % linear array test %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Find P for a linear array
                phi_c = obj.Speaker_Array_Centre/180*pi;
                C = O + obj.res*Rl*[cos(phi_c) -sin(phi_c) 0;...
                    sin(phi_c)  cos(phi_c) 0;...
                    0           0      1]*[1;0;0];
                m1 = [cos(theta);sin(theta);0];
                m2 = [cos(phi_c-pi/2);sin(phi_c-pi/2);0];
                A=[m1 -m2];
                b = [C-B];
                w = pinv(A)*b;
                P = B+w(1)*m1;
                P_ = C+w(2)*m2;
                if ~all(round(P,10)==round(P_,10))
                    error('Planewave does not intersect loudspeaker array.');
                end
                
                Pv = P-O;
                PQv = Qv-Pv;
                d = sqrt(sum((Q-P).^2))/obj.res;
                Theta = abs(pi - abs(theta) - abs(phi_c) );
                d_alias=[];
                gam_found=[];
                for pm = [-1 1]
                    % maximum allowable angle of grating lobe
                    alph_=pm*asin((r1+r2)/d);
                    PAv=[cos(alph_) -sin(alph_) 0;...
                        sin(alph_)  cos(alph_) 0;...
                        0           0      1] * PQv;
                    
                    hold on; plot3([P(1) PAv(1)+P(1)],...
                        [P(2) PAv(2)+P(2)],...
                        [P(3) PAv(3)+P(3)],'--k'); hold off;
                    
                    d_alias(end+1) = abs(  Pv.'   * (rotMat_90*PAv  )) ...
                        / sqrt(sum(( PAv ).^2)) / obj.res;
                    
                    
                    % grating lobe angle
                    %                     gam_ = pm * asin(2*pi/obj.k_global/(obj.Speaker_Array_Length/(L-1)));
                    PBv = Bv-Pv;
                    
                    %tmp = asin(r1 / (sqrt(sum(( PBv ).^2)) / obj.res));
                    %                     tmp=0;
                    gam_ = pm * ( asin(  2*pi*(L-1) / (k*arrLen) - sin(Theta) )  + Theta);
                    
                    gv=[cos(gam_) -sin(gam_) 0;...
                        sin(gam_)  cos(gam_) 0;...
                        0           0      1] * PBv;
                    
                    hold on; plot3([P(1) gv(1)+P(1)],...
                        [P(2) gv(2)+P(2)],...
                        [P(3) gv(3)+P(3)],':k'); hold off;
                    
                    %                     d_g(end+1) = abs(  Pv.'   * (rotMat_90*PAv  )) ...
                    %                         / sqrt(sum((PAv  ).^2)) / obj.res;
                    gam_found(end+1) = atan2(norm(cross(PAv,PBv)), dot(PAv,PBv)); % EQUIVALENT MATHEMATICALLY --> acos(dot(v1, v2) / (norm(v1) * norm(v2)))
                    
                end
                
                d_alias = max(d_alias);
                gam_found = max(gam_found);
                Alias_k =   2*pi*(L-1) / ((sin(gam_found - Theta) + sin(Theta))*arrLen)  ;
                disp(['Alias Frequency: ' num2str(round(Alias_k/2/pi*obj.c,0)) 'Hz']);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                plotCircle(obj, zeros(1,3), P/obj.res, ...
                    obj.Multizone_Soundfield.Bright_Zone.Radius_q, ...
                    {'--', 'Color', [0 0 0], 'LineWidth', LineWid});
                
                
                
            elseif strcmpi('circle',obj.Speaker_Array_Type)
                
                alpha_ = pi +theta - asin(Br*sin(-theta + Bth)/Rl);
                % Find tangent
                P = O + Rl*obj.res * ...
                    [cos(alpha_) -sin(alpha_) 0;...
                    sin(alpha_)  cos(alpha_) 0;...
                    0           0       1]*[1;0;0];
                Pv = P-O;
                PQv = Qv-Pv;
                d = sqrt(sum((Q-P).^2))/obj.res;
                d_alias=[];
                for pm = [-1 1]
                    
                    alph_=pm*asin((r1+r2)/d);
                    
                    PAv=[cos(alph_) -sin(alph_) 0;...
                        sin(alph_)  cos(alph_) 0;...
                        0           0      1] * PQv;
                    hold on; plot3([P(1) PAv(1)+P(1)],...
                        [P(2) PAv(2)+P(2)],...
                        [P(3) PAv(3)+P(3)],'--k'); hold off;
                    
                    %                 P_r1 = rotMat_90*PAv / sqrt(sum(PAv.^2)) * r1*obj.res;
                    %                 P_1 = P + P_r1;
                    %                 P_2 = P_1 + PAv;
                    %                 hold on; plot3([P_1(1) P_2(1)],...
                    %                     [P_1(2) P_2(2)],...
                    %                     [P_1(3) P_2(3)]*0,'--k'); hold off;
                    
                    %             d_found  = sum((O - P_1) .* P_r1) / sqrt(sum(P_r1.^2)) / obj.res;
                    %             d_found2 = sum((O - P)   .* P_r1) / sqrt(sum(P_r1.^2)) / obj.res;
                    d_alias(end+1) = abs(  Pv.'   * (rotMat_90*PAv  )) ...
                        / sqrt(sum((PAv  ).^2)) / obj.res;
                end
                d_alias = max(d_alias);
                
                plotCircle(obj, zeros(1,3), O/obj.res, d_alias, ...
                    {'--', 'Color', [1 0 0], 'LineWidth', LineWid});
                
                
                % Perpendicular distance from planewave to center
                PBv = Bv-Pv;
                d_pb = abs(  Pv.'   * (rotMat_90*PBv  )) ...
                    / sqrt(sum((PBv  ).^2)) / obj.res;
                
                % Maximum allowable aliasing frequency
                Alias_k = max( ...
                    [(2*pi*L - phi_L) / ((d_alias + d_pb)*phi_L), ...
                    (2*pi*L - phi_L) / (2*Rprime*phi_L)]);
                disp(['Alias Frequency: ' num2str(round(Alias_k/2/pi*obj.c,0)) 'Hz']);
                
                % Current grating lobe edge
                Alias_R = (2*pi*L - phi_L) / (obj.k_global*phi_L) - d_pb;
                plotCircle(obj, zeros(1,3), O/obj.res, Alias_R, ...
                    {'--', 'Color', [1 0 1], 'LineWidth', LineWid});
                
                %             M =  tan(pi/2 + asin( Rl*sin(Qth-alpha_) / d) - asin((r1+r2)/d)) ; % Gradient of aliasing tangent
                %             x = 0:300;
                %             % Maximum allowable aliasing lobe center
                %             y = M*(x-Alias_Point(1))+Alias_Point(2);
                %             hold on; plot3(x,y,x*0+maxZ,'-.','Color', [0 0 0], 'LineWidth', LineWid); hold off;
                %             % Maximum allowable aliasing lobe edge
                %             y = M*(x-Alias_Point(1) + sin(atan(M))*r1*obj.res)+Alias_Point(2) + cos(atan(M))*r1*obj.res;
                %             hold on; plot3(x,y,x*0+maxZ,':','Color', [0 0 0], 'LineWidth', LineWid); hold off;
                %
                %             % Maximum allowable aliasing radius
                %             Alias_R_found = abs( O(2) - M*O(1) ...
                %                 + M*Alias_Point(1) - M*sin(atan(M))*r1*obj.res - Alias_Point(2) - cos(atan(M))*r1*obj.res ) ...
                %                 / sqrt(1^2 + M^2) / obj.res;
            end
            
            plotCircle(obj, zeros(1,3), P/obj.res, ...
                obj.Multizone_Soundfield.Bright_Zone.Radius_q, ...
                {'--', 'Color', [0 0 0], 'LineWidth', LineWid});
            
            % Virtual point
            arrow(P,BZ_Centre,arrowSpecs{:});
            
            
            
            
            
            
            
            
            
            
            
            
            hold off
        end
        
        function drawPlaneWaveArrow(obj,startPt,endPt,PWangle,arrowSpecs)
            P = endPt-startPt;
            [~,el,r]=cart2sph(P(1),P(2),P(3));
            r_arr = r*2;
            [P(1),P(2),P(3)]=sph2cart(PWangle/180*pi,el,r_arr);
            [Pmid(1),Pmid(2),Pmid(3)]=sph2cart(PWangle/180*pi,el,r_arr/2);
            endPt = P + startPt;
            arrow(startPt,endPt,arrowSpecs{:});
            %draw three perpendicular straight lines
            noLines = 4;
            lineLength = 7;
            sepDist = 10;
            x = repmat(linspace(-1,1,noLines),3,1)*sepDist;
            y = repmat([-1;0;1],1,noLines)*lineLength;
            [az,~,r]=cart2sph(x,y,zeros(size(x)));
            [Pts(:,:,1),Pts(:,:,2),Pts(:,:,3)] = sph2cart(az+PWangle/180*pi,ones(size(az))*el,r);
            Pts = Pts + permute(repmat(Pmid'+startPt,1,noLines,3),[3 2 1]);
            hold on;
            plot3(gca, Pts(:,:,1),Pts(:,:,2),Pts(:,:,3), 'Color', [0 0 0], 'LineWidth', 2);
            hold off;
        end
        
        function addText2Arrow(obj,txt,startPt,endPt,fontsz,DistAway,AngAway)
            P=endPt-startPt;
            [az,el,r] = cart2sph(P(1),P(2),P(3));
            [x,y,z] = sph2cart(az,el,r/2);
            [x2,y2,z2] = sph2cart(az+pi/2-AngAway/180*pi,0,DistAway);
            P=[x;y;z]+[x2;y2;z2]+startPt;
            
            text(P(1),P(2),P(3), ...
                ['\fontsize{' num2str(fontsz) '} \fontname{times}' txt], ...
                'interpreter','tex', ...
                'HorizontalAlignment','center')
        end
        
        function drawAngle(obj,centre,startAng,endAng,Radius,LineWid,Elevation)
            az = startAng : ((endAng>startAng)*2-1)*pi/obj.res : endAng;
            [x,y,z]=sph2cart(az,zeros(size(az))*Elevation,Radius);
            P = [x;y;z] + repmat(centre,1,length(az));
            hold on;
            plot3(gca, P(1,:), P(2,:), P(3,:), 'Color', [0 0 0], 'LineWidth', LineWid);
            hold off;
        end
        
        function addText2Angle(obj,txt,centre,startAng,endAng,fontsz,DistAway,Elevation)
            az = median( startAng : ((endAng>startAng)*2-1)*pi/obj.res : endAng );
            [x,y,z]=sph2cart(az,Elevation,DistAway);
            P = [x;y;z] + centre;
            
            text(P(1),P(2),P(3), ...
                ['\fontsize{' num2str(fontsz) '} \fontname{times}' txt], ...
                'interpreter','tex', ...
                'HorizontalAlignment','center')
        end
        
        function h = sourcePlots(obj, field, realistic, details)
            if nargin < 4; details = []; end
            if nargin < 3; realistic = false; end
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            if realistic
                maxZ = 0;
            else
                maxZ = max(real(field(:)));
            end
            %O = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
            O = [flip(size(field))/2 0] + [flip(obj.Origin)*obj.res 0];
            
            hold on;
            
            [sl_x, sl_y] = pol2cart(obj.Loudspeaker_Locations(:,1), (obj.Loudspeaker_Locations(:,2) * obj.res));
            
            if strcmpi( obj.Multizone_Soundfield.Bright_Zone.SourceType, 'ps')
                ps_d = obj.Multizone_Soundfield.Bright_Zone.SourceOrigin.Distance;
                ps_a = obj.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle/180*pi;
                [ps_x, ps_y] = pol2cart(ps_a, (ps_d * obj.res)-1);
            else
                ps_x = nan;
                ps_y = nan;
            end
            
            if realistic
                L = obj.Loudspeaker_Count;
                widL = repmat(obj.Loudspeaker_Dimensions(1) * obj.res, L, 1 ); %metres
                depL = repmat(obj.Loudspeaker_Dimensions(2) * obj.res, L, 1 ); %metres
                higL = repmat(obj.Loudspeaker_Dimensions(3) * obj.res, L, 1 ); %metres
                widL_half = widL / 2; %metres
                depL_half = depL / 2; %metres
                higL_half = higL / 2; %metres
                
                thL = zeros(L,4);
                xL = thL;
                yL = thL;
                radL = zeros(L,1);
                
                if strcmp(obj.Speaker_Array_Type, 'circle')
                    thL = repmat(obj.Loudspeaker_Locations(:,1),1,4);
                    radL = repmat(obj.Radius * obj.res,L,1) - 1;
                    [xL, yL] = pol2cart(obj.Loudspeaker_Locations(:,1), ...
                        obj.Loudspeaker_Locations(:,2) * obj.res - 1);
                    xL = repmat(xL,1,4);
                    yL = repmat(yL,1,4);
                elseif ~isempty(strfind(obj.Speaker_Array_Type, 'line'))
                    thL = repmat(obj.Speaker_Array_Centre,L,4)/180*pi;
                    [xL, yL] = pol2cart(obj.Loudspeaker_Locations(:,1), ...
                        obj.Loudspeaker_Locations(:,2) * obj.res - 1);
                    xL = repmat(xL,1,4);
                    yL = repmat(yL,1,4);
                end
                
                
                for i = 1:obj.Loudspeaker_Count
                    % Create Cube
                    px_=[-1,1,1,-1;-1,1,1,-1];
                    py_=[-1,-1,1,1;-1,-1,1,1];
                    pz_=[-1,-1,-1,-1;1,1,1,1];
                    %Size to Rect Prism
                    px = ([px_;py_;pz_]'+1) .* widL_half(i); %plus one shifts cube so origin is on a face
                    py = ([py_;pz_;px_]'+0) .* depL_half(i);
                    pz = ([pz_;px_;py_]'+0) .* higL_half(i);
                    %Rotate Rect Prism
                    [az, el, r] = cart2sph(px, py, pz);
                    [px, py, pz] = sph2cart(az + thL(i,1), el + 0, r + 0);
                    %Shift Rect Prism
                    px = px + O(1) + xL(i,1);
                    py = py + O(2) + yL(i,1);
                    %pz = pz + O + zL(i,1);
                    
                    h.box(i) = patch(px,py,pz,[0.5 0.5 0.5]);
                    
                end
            end
            [slx,sly,slz] = sphere(5);
            if isempty(details)
                sphereRad = 2; %centimetres
            else
                sphereRad = details.sphereRad;
            end
            sphereRad = sphereRad/100*obj.res;
            for i = 1:obj.Loudspeaker_Count
                h.sphere(i) = surf( ...
                    slx*sphereRad + sl_x(i) + O(1), ...
                    sly*sphereRad + sl_y(i) + O(2), ...
                    slz*sphereRad + ones(size(sl_x(i)))*maxZ, ...
                    repmat(0.5,[size(slx),3]), ...
                    'linestyle',':');
            end
            %scatter3(sl_x + O, sl_y + O, ones(length(sl_x),1)*maxZ, 100, 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', [1 1 1], 'LineWidth', 3);
            
            % Plot virtual source
            if strcmpi( obj.Multizone_Soundfield.Bright_Zone.SourceType, 'ps')
                BZ_O =[obj.Multizone_Soundfield.Bright_Zone.Origin_q.X, ...
                    obj.Multizone_Soundfield.Bright_Zone.Origin_q.Y]*obj.res;
                h.sphereVSrc = surf( ...
                    slx*sphereRad + ps_x + O(1) + BZ_O(1), ...
                    sly*sphereRad + ps_y + O(2) + BZ_O(2), ...
                    slz*sphereRad + ones(size(ps_x))*maxZ, ...
                    repmat(0.5,[size(slx),3]), ...
                    'linestyle','none');
                h.sphereVSrc.CData = permute(repmat([0 0.2 1].',[1,size(slx)]),[2 3 1]);
            end
            [c_x, c_y] = pol2cart(obj.Speaker_Array_Centre/180*pi, (obj.Radius * obj.res));
            
            h.sphereVSrc = surf( ...
                slx*sphereRad + c_x + O(1), ...
                sly*sphereRad + c_y + O(2), ...
                slz*sphereRad + ones(size(ps_x))*maxZ, ...
                repmat(0.5,[size(slx),3]), ...
                'linestyle','none');
            h.sphereVSrc.CData = permute(repmat([1 0 0].',[1,size(slx)]),[2 3 1]);
            
            
            hold off;
        end
        
    end
    
end

