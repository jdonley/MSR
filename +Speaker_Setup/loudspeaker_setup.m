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
        res = 50;  % samples per metre %Resolution of soundfield
        Loudspeaker_Count = 32;     % Number of speakers used in the reproduction (Number of speakers required = ceil( Angular_Window * (2*M + 1) / (2*pi) ) ).
        Radius = 1.5;               % Radius of speaker layout
        Speaker_Arc_Angle = 180;    % Angle of the arc of speakers (360 = full circle, 180 = semi-circle)
        Angle_FirstSpeaker = 0;     % The angle where the first speaker of the array occurs
        Speaker_Array_Centre = 180;   % Angle of the centre of the loudspeaker arc
        Speaker_Spacing = 0.01;       % The spacing between each consecutive loudspeaker
        Speaker_Array_Type = 'circle'; % 'circle' or 'line'
        Loudspeaker_Dimensions;
        Loudspeaker_Type = 'Genelec 8010A';
        k_global = 2000 /343*2*pi;  % Frequency in wavenumber
        Multizone_Soundfield;       % multizone_soundfield object for reproduction
        %N_zone_samples = 1;             % The number of samples in each zone
        %Sampling_method = 'Centre';            % The layout of the samples in the zones. ('Centre', 'Random', 'Square', 'Circlular', 'All')
        
        % Results
        % 3D
        Loudspeaker_Weights = [];
        Loudspeaker_Locations = [];
        Soundfield_reproduced = [];    % The complex values of the reproduced sound field that is produced from this class.
        Desired_Mask = [];
        %         Soundfield_Error = [];         % The error in the bright field
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
    end
    
    properties (Access = private)
        
    end
    
    %% Contructor
    methods
        function obj = multizone_soundfield_OBE(~)
            
        end
    end
    
    %% Private Methods
    methods (Access = private)
        
        function obj = createEmptySoundfield(obj)
            width = int16(obj.res * obj.Radius * 2);
            obj.Soundfield_reproduced = zeros(width,width);
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
                elseif nargin < 3
                    obj = obj.setRadius(obj.Radius);
                else
                    obj = obj.setRadius(Radius);
                end
                
            obj = obj.setWavenumberFromChild();
            
            obj = obj.createEmptySoundfield;            
                       
            obj = obj.save_Bright_Samples();
            obj = obj.save_Quiet_Samples();
            
            if ~strcmp(Debug, 'SAMPLES_ONLY')
                O = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
                th1 = 0;
                r1 = 0;
                th2 = obj.Loudspeaker_Locations(:,1) + pi;
                r2 = obj.Loudspeaker_Locations(:,2) * obj.res;
                H = 0;
                L = obj.Loudspeaker_Weights;
                                
                xy = O - (1:O*2);
                wid = length(obj.Soundfield_reproduced);
                Nspkrs = obj.Loudspeaker_Count;
                
                [xx,yy] = ndgrid(xy);
                th1 = atan2(xx, yy);
                r1 = sqrt( xx.^2 + yy.^2 );
                
                % Create matrices
                r1_  = repmat(r1 ,1,1,Nspkrs);
                th1_ = repmat(th1,1,1,Nspkrs);
                r2_  = permute( repmat(r2     ,1,wid,wid), [2 3 1] );
                th2_ = permute( repmat(th2    ,1,wid,wid), [2 3 1] );
                L_   = permute( repmat(conj(L),1,wid,wid), [2 3 1] );
                
                r = (r1_.^2 + r2_.^2 - 2*r1_.*r2_.*cos(th1_-th2_)).^(1/2) / obj.res;
                
                H = permute( 1i/4 * besselh(0, obj.k_global * r ), [2 1 3]);

                obj.Soundfield_reproduced = sum( L_.*H , 3 )';
                
                %obj = obj.norm_soundfield();
            end
            
            % Equal delay 
            obj.Loudspeaker_Weights     = obj.Loudspeaker_Weights * abs( besselh(0, obj.k_global ) ); %Remove effect of frequency dependent decay
            % Normalise to the maximum of the absolute bright zone values
            maxBrightVal = mean(abs(obj.Bright_Samples(:))); % Changed to find the average as our goal is to have the bright zone fit on average.
            obj.Loudspeaker_Weights     = obj.Loudspeaker_Weights   / maxBrightVal; %Added to keep speaker weights tied with their response
            obj.Bright_Samples          = obj.Bright_Samples        / maxBrightVal;
            obj.Quiet_Samples           = obj.Quiet_Samples         / maxBrightVal;
            obj.Bright_Sample           = obj.Bright_Sample         / maxBrightVal;
            obj.Quiet_Sample            = obj.Quiet_Sample          / maxBrightVal;
            obj.Soundfield_reproduced   = obj.Soundfield_reproduced / maxBrightVal;
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
            
            Q0 = obj.Loudspeaker_Count; 
            phi = obj.Speaker_Arc_Angle / 180 * pi;
            delta_phi_s = phi / Q0;
            
            phi_q = obj.Loudspeaker_Locations(:,1) + pi;
            R_q = obj.Loudspeaker_Locations(:,2);            
            M = obj.Multizone_Soundfield.getGlobalModeLimit;            
            
            Alpha_to_Beta_Coeff = zeros(Q0,2*M+1);
            for m = -M:M
                Alpha_to_Beta_Coeff(:,m+M+1) = 2 ./ (1i*pi*besselh(m, obj.Multizone_Soundfield.k_global * R_q));
            end
            Beta =   Alpha_to_Beta_Coeff .* repmat(obj.Multizone_Soundfield.Alpha_Coeffs, Q0, 1) ;
            
            m = -M:M;
            obj.Loudspeaker_Weights = zeros(1, Q0);
            
            for q = 1:Q0
                obj.Loudspeaker_Weights(q) = sum(Beta(q,:) .* exp(1i*m*phi_q(q,:)) * delta_phi_s);
            end
            obj.Loudspeaker_Weights = obj.Loudspeaker_Weights';
        end
        
        
        function obj = calc_Loudspeaker_Locations(obj)        
            Q0 = obj.Loudspeaker_Count;
            Rl = obj.Radius;
            
            if strcmp(obj.Speaker_Array_Type, 'circle')
                phi = obj.Speaker_Arc_Angle / 180 * pi;
                delta_phi_s = phi / Q0;
                first_spkr_offset = obj.Angle_FirstSpeaker / 180 * pi;
                phi_q = ((1:Q0) - 1) * delta_phi_s + pi + first_spkr_offset;
                Rl_ = repmat(Rl, length(phi_q), 1);
                obj.Loudspeaker_Locations = [(phi_q(:) - pi) Rl_];
                
            elseif strcmp(obj.Speaker_Array_Type, 'line')
                wid = obj.Loudspeaker_Dimensions(1);
                space = obj.Speaker_Spacing;
                centre = obj.Speaker_Array_Centre;
                L = obj.Loudspeaker_Count;
                len = (wid+space)*L - space;
                x = repmat(obj.Radius, 1, L);
                
                if iseven(L)
                    y = ((wid+space)/2):(wid+space):len/2;
                    y = [flip(-y), y];
                else
                    y = (wid+space):(wid+space):len/2 + (wid+space)/2;
                    y = [flip(-y), 0, y];
                end
                
                [spkr_theta, spkr_rho]=cart2pol(x,y);
                spkr_theta = spkr_theta + centre/180*pi;
                
                obj.Loudspeaker_Locations = [spkr_theta(:), spkr_rho(:)];
                
            elseif strcmp(obj.Speaker_Array_Type, 'coprime')
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
        
        function obj = calcDesiredMask(obj)
            width = int16(obj.res * obj.Radius * 2);
            xyvec = linspace(-obj.Radius, obj.Radius , width);
            [xx,yy] = meshgrid(xyvec,xyvec);
            obj.Desired_Mask = xx.^2 + yy.^2 <= (obj.Multizone_Soundfield.Radius)^2 * ones(width, width);
        end
        
        function f = getFrequency(obj)
            f = obj.k_global * 343 / (2 * pi);
        end
        
        function obj = setWavenumberFromChild(obj)
           obj.k_global = obj.Multizone_Soundfield.k_global;
        end
        
        function obj = setRadius(obj, Radius)
           obj.Radius = Radius;           
           obj = obj.calcDesiredMask();
        end
                        
        
        function obj = setLoudspeakerType(obj, Type)
            obj.Loudspeaker_Type = Type;
            
            if strcmp(Type, 'Genelec 8010A')
                obj.Loudspeaker_Dimensions = [0.121, 0.116, 0.195]; %Genelec 8010A
                
            elseif strcmp(Type, 'Genelec 8020C')
                obj.Loudspeaker_Dimensions = [0.151, 0.142, 0.242]; %Genelec 8020C
                
            elseif strcmp(Type, 'Meyer MM-4XP')
                obj.Loudspeaker_Dimensions = [0.103, 0.145, 0.103]; %Meyer MM-4XP
                
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
                warning(['The previously specified loudspeaker arc is too small to accommodate the desired number of loudspeakers.' 10, ...
                         'Either:' 10 ...
                         9	'1) Reduce the number of loudspeakers.' 10 ...
                         9  '2) Use smaller loudspeakers.  or' 10 ... 
                         9  '3) Increase the size of the loudspeaker arc.' 10, ...
                         'The loudspeaker arc will be adjusted now...' 10, ...
                         9  'Old loudspeaker arc: ~' num2str(round(obj.Speaker_Arc_Angle,1)) '°' 10, ...
                         9  'New loudspeaker arc: ~' num2str(round(new_Speaker_Arc_Angle,1)) '°' 10 10]);
            end
            obj.Speaker_Arc_Angle = new_Speaker_Arc_Angle;
            obj.Angle_FirstSpeaker = centre - obj.Speaker_Arc_Angle/2 + (((space + width)/2) / obj.Radius * 180/pi );
        end
        
        
        function obj = norm_soundfield(obj)
            f_ = obj.Soundfield_reproduced .* obj.Desired_Mask;
            
            obj.Soundfield_reproduced = obj.Soundfield_reproduced / max(abs(real(f_(:)))) * 4/3;
        end
        
        function obj = save_Bright_Samples(obj)
                        
%             % If number is not a square number change method to random
%             if ~~rem(obj.N_zone_samples, sqrt(obj.N_zone_samples))
%                 obj.Sampling_method = 'Random';
%             else
%                 
%             end
%             
%             if strcmp(obj.Sampling_method, 'Random')
%                 
%             end

O = obj.Multizone_Soundfield.Bright_Zone.Radius_q * obj.res;
Oz = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
th2 = obj.Loudspeaker_Locations(:,1) + pi;
r2 = obj.Loudspeaker_Locations(:,2) * obj.res;
L = obj.Loudspeaker_Weights;

obj.Bright_Samples_Locations = ones(O*2, O*2, 2)*NaN;           
            
            xy = O - (1:O*2);
            wid = O*2;
            Nspkrs = obj.Loudspeaker_Count;
            mask = obj.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask;
            x_off = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Bright_Zone.Origin_q.Y * obj.res;
            y_off = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Bright_Zone.Origin_q.X * obj.res;
            
            [xx,yy] = ndgrid(xy);
            xx = Oz - (xx .* mask + x_off);
            yy = Oz - (yy .* mask + y_off);
            th1 = atan2(xx, yy);
            r1 = sqrt( xx.^2 + yy.^2 );
            obj.Bright_Samples_Locations(:,:,1) = th1 + pi;
            obj.Bright_Samples_Locations(:,:,2) = r1 / obj.res;
            
            
            % Create matrices
            r1_  = repmat(r1 ,1,1,Nspkrs);
            th1_ = repmat(th1,1,1,Nspkrs);
            r2_  = permute( repmat(r2     ,1,wid,wid), [2 3 1] );
            th2_ = permute( repmat(th2    ,1,wid,wid), [2 3 1] );
            L_   = permute( repmat(conj(L),1,wid,wid), [2 3 1] );
            
            r = (r1_.^2 + r2_.^2 - 2*r1_.*r2_.*cos(th1_-th2_)).^(1/2) / obj.res;
            
            H = permute( 1i/4 * besselh(0, obj.k_global * r ), [2 1 3]);
            
            obj.Bright_Samples = sum( L_.*H , 3 );
            
            obj.Bright_Sample = mean( abs( obj.Bright_Samples(:) ) );
        end
        
        function obj = save_Quiet_Samples(obj)
            
%             % If number is not a square number change method to random
%             if ~~rem(obj.N_zone_samples, sqrt(obj.N_zone_samples))
%                 obj.Sampling_method = 'Random';
%             else
%                 
%             end
%             
%             if strcmp(obj.Sampling_method, 'Random')
%                 
%             end

O = obj.Multizone_Soundfield.Quiet_Zone.Radius_q * obj.res;
Oz = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
th2 = obj.Loudspeaker_Locations(:,1) + pi;
r2 = obj.Loudspeaker_Locations(:,2) * obj.res;
L = obj.Loudspeaker_Weights;

obj.Quiet_Samples_Locations = ones(O*2, O*2, 2)*NaN;


            xy = O - (1:O*2);
            wid = O*2;
            Nspkrs = obj.Loudspeaker_Count;
            mask = obj.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask;
            x_off = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Quiet_Zone.Origin_q.Y * obj.res;
            y_off = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Quiet_Zone.Origin_q.X * obj.res;
            
            [xx,yy] = ndgrid(xy);
            xx = Oz - (xx .* mask + x_off);
            yy = Oz - (yy .* mask + y_off);
            th1 = atan2(xx, yy);
            r1 = sqrt( xx.^2 + yy.^2 );
            obj.Quiet_Samples_Locations(:,:,1) = th1 + pi;
            obj.Quiet_Samples_Locations(:,:,2) = r1 / obj.res;
            
            
            % Create matrices
            r1_  = repmat(r1 ,1,1,Nspkrs);
            th1_ = repmat(th1,1,1,Nspkrs);
            r2_  = permute( repmat(r2     ,1,wid,wid), [2 3 1] );
            th2_ = permute( repmat(th2    ,1,wid,wid), [2 3 1] );
            L_   = permute( repmat(conj(L),1,wid,wid), [2 3 1] );
            
            r = (r1_.^2 + r2_.^2 - 2*r1_.*r2_.*cos(th1_-th2_)).^(1/2) / obj.res;
            
            H = permute( 1i/4 * besselh(0, obj.k_global * r ), [2 1 3]);
            
            obj.Quiet_Samples = sum( L_.*H , 3 );
                        
            obj.Quiet_Sample = mean( abs( obj.Quiet_Samples(:) ) );
        end
        
        %% Plotting Functions
        function h = plotSoundfield(obj, field, colour, realistic)
            if nargin < 4; realistic = false; end;
            if nargin < 3; colour = 'default'; end;
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            ZaxisSize = max(abs(real(field(:))));
            f = real(field) .* obj.Desired_Mask;
            CaxisSize = max(abs(f(:))) * 3/4;
            
            XYTick = [1 length(field)/4 length(field)/2 length(field)*3/4 length(field)];
            XYTickLabel = [ -length(field)/2/obj.res, -length(field)/4/obj.res, 0, length(field)/4/obj.res, length(field)/2/obj.res];
            
            
            h = surf(real(field),'EdgeColor','None');
            if realistic
                drawnow; pause(0.05);
                h.FaceAlpha = 0.5;
            end
            colormap(colour);
            view(2);
            title('Real',...
                'FontWeight','bold',...
                'FontSize',24,...
                'FontName','Arial');
            axis equal;
            box on;
            if realistic
                room = ([4.815, 3.3]-0.10)*obj.res;  %The subtraction is the wall thickness compensation for the anechoic chamber
                axis([(size(field,1)/2-room(1)/2) (size(field,1)/2+room(1)/2) (size(field,2)/2-room(2)/2) (size(field,2)/2+room(2)/2)]);
            else
                axis([0 size(field,1)+1 0 size(field,2)+1]);
            end
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
            
            obj.zonePlots(real(field),realistic);
            obj.speakerPlots(real(field),realistic);
            
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
        
        function zonePlots(obj, field, realistic)
            if nargin < 3; realistic = false; end
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            % For a Global Zone outline plot
            GlobalZone_temp = obj.Multizone_Soundfield.Bright_Zone;
            GlobalZone_temp.Origin_q = struct('X',0,'Y',0);
            GlobalZone_temp.Radius_q = obj.Multizone_Soundfield.Radius;
            
            if realistic
                maxZ = 0;
            else
                maxZ = max(real(field(:)));
            end
            hold on
            th = 0:pi/obj.res:2*pi;
            SpatialZones = [obj.Multizone_Soundfield.Quiet_Zone obj.Multizone_Soundfield.Bright_Zone GlobalZone_temp];
            
            for q = 1:length(SpatialZones)
                Ox = SpatialZones(q).Origin_q.X * obj.res + length(field) / 2 + 1;
                Oy = SpatialZones(q).Origin_q.Y * obj.res + length(field) / 2 + 1;
                xunit = SpatialZones(q).Radius_q * obj.res * cos(th) + Ox;
                yunit = SpatialZones(q).Radius_q * obj.res * sin(th) + Oy;
                plot3(xunit, yunit, ones(length(xunit),1)*maxZ, 'Color', [0 0 0], 'LineWidth', 3);
            end            
            
            hold off
        end
        
        function speakerPlots(obj, field, realistic)
            if nargin < 3; realistic = false; end
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            if realistic
                maxZ = 0;
            else
                maxZ = max(real(field(:)));
            end
            O = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
            hold on;            
            
            [sl_x, sl_y] = pol2cart(obj.Loudspeaker_Locations(:,1), (obj.Loudspeaker_Locations(:,2) * obj.res)-1);
            
            if realistic                
                L = obj.Loudspeaker_Count;                
                widL = repmat(obj.Loudspeaker_Dimensions(1) * obj.res, L, 1 ); %metres
                depL = repmat(obj.Loudspeaker_Dimensions(2) * obj.res, L, 1 ); %metres
                widL_half = widL / 2; %metres
                
                thL = zeros(L,4);
                xL = thL;
                yL = thL;
                radL = zeros(L,1);
                
                if strcmp(obj.Speaker_Array_Type, 'circle')
                    thL = repmat(obj.Loudspeaker_Locations(:,1),1,4);
                    radL = repmat(obj.Radius * obj.res,L,1) - 1;
                elseif strcmp(obj.Speaker_Array_Type, 'line')
                    thL = repmat(obj.Speaker_Array_Centre,L,4)/180*pi;
                    [xL, yL] = pol2cart(obj.Loudspeaker_Locations(:,1), ...
                        obj.Loudspeaker_Locations(:,2) * obj.res - 1);
                    xL = repmat(xL,1,4);
                    yL = repmat(yL,1,4);
                end
                
                X = [radL, radL, radL + depL, radL + depL ];
                Y = [widL_half, -widL_half, -widL_half, widL_half ];                
                [TH, R] = cart2pol(X, Y);                
                [X, Y] = pol2cart(TH + thL, R);                
                X = X + xL;
                Y = Y + yL;
                
                for i = 1:obj.Loudspeaker_Count
                    spkr_crnrs_x = [X(i,:), X(i,1)] + O;
                    spkr_crnrs_y = [Y(i,:), Y(i,1)] + O;
                    spkr_Z = repmat(maxZ, 1, 5);
                    plot3(spkr_crnrs_x, spkr_crnrs_y, spkr_Z, 'Color', [0 0 0], 'LineWidth', 2);
                end
            end
            scatter3(sl_x + O, sl_y + O, ones(length(sl_x),1)*maxZ, 100, 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', [1 1 1], 'LineWidth', 3);
            
            
            hold off;
        end
        
    end
    
end

