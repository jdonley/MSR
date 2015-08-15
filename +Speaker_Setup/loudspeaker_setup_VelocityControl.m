classdef loudspeaker_setup_VelocityControl
    % LOUDSPEAKER_SETUP_VELOCITY_CONTROL - This class provides practical setup information
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
        k_global = 2000 /343*2*pi;  % Frequency in wavenumber
        Multizone_Soundfield;       % multizone_soundfield object for reproduction
        %N_zone_samples = 1;             % The number of samples in each zone
        %Sampling_method = 'Centre';            % The layout of the samples in the zones. ('Centre', 'Random', 'Square', 'Circlular', 'All')
        
        % Results
        % 3D
        Loudspeaker_Weights_Pressure = [];
        Loudspeaker_Weights_Velocity = [];
        Loudspeaker_Locations = [];
        Soundfield_reproduced = [];    % The complex values of the reproduced sound field that is produced from this class.
        Soundfield_PressureReproduced = [];
        Soundfield_VelocityReproduced = [];
        Intensity_Soundfield = [];
        
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
        
        transferFunc_Pressure = [];
        transferFunc_Velocity = [];
        
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
%             else
%                 if nargin < 2
%                     obj = obj.setRadius(obj.Radius);
%                     Debug = '';
%                 elseif nargin < 3
%                     obj = obj.setRadius(obj.Radius);
%                     obj.N_zone_samples = 1;
%                 else
%                     obj.N_zone_samples = Radius;
%                 end
%             end
            
            obj = obj.calc_Loudspeaker_Weights(Debug);
            
            obj = obj.createEmptySoundfield;            
                                  
            if ~strcmp(Debug, 'SAMPLES_ONLY')
                O = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
                th1 = 0;
                r1 = 0;
                P1 = obj.Loudspeaker_Locations;                
                H = 0;
                Z = 0;
                Lp = obj.Loudspeaker_Weights_Pressure;
                Lv = obj.Loudspeaker_Weights_Velocity;
                
                n=0;
                for x = 1:O*2
                    for y = 1:O*2
                        P2 = repmat([x-O, y-O],[obj.Loudspeaker_Count 1]) / obj.res;
                        
                        d =  sqrt( sum( (P2 - P1).^2, 2 ) );
                                                
                        k = obj.k_global;
                        omega = 343 * k;
                        rho_0 = 20*1e-6;
                        
                        Zp = 1j * omega * rho_0 * exp(-1j*obj.k_global*d) / 4 / pi ./ d;
                        obj.Soundfield_PressureReproduced(y, x) =  Lp' * Zp; % Sum pressures
                        
                        Zv = 1j*k*exp(-1j*k*d)./(4*pi*d) .* (1 + 1./(1j*k*d));
                        obj.Soundfield_VelocityReproduced(y, x) =  Lv' * Zv;
                    end
                    if strcmp('DEBUG', Debug)
                        fprintf(repmat('\b',1,n));
                        n=fprintf('%.2f%%', x / (O*2) * 100);
                    end
                end
                
                %obj = obj.norm_soundfield();
            end
            obj.Soundfield_reproduced = obj.Soundfield_VelocityReproduced;
            obj.Soundfield_reproduced(isinf(obj.Soundfield_reproduced))=0;
            %surf(real(obj.Soundfield_reproduced),'LineStyle','none');view(2);
            
            % Normalise to the maximum of the absolute bright zone values
            obj = obj.save_Bright_Samples();
            obj = obj.save_Quiet_Samples();
            maxBrightVal = max(abs(obj.Bright_Samples(:)));
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
            
            phi = obj.Speaker_Arc_Angle / 180 * pi;
            first_spkr_offset = obj.Angle_FirstSpeaker / 180 * pi;
           
            %obj.Loudspeaker_Count = ceil( phi_window*(2*M+1)/(2*pi) );
            Q0 = obj.Loudspeaker_Count; 
            q = 1:Q0;
            Rl = repmat(obj.Radius,size(q)); 
                       
             delta_phi_s = phi / Q0;
             phi_q = ((1:Q0) - 1) * delta_phi_s + pi + first_spkr_offset;
         
             [x, y] = pol2cart(phi_q - pi, Rl);
             obj.Loudspeaker_Locations = [x', y']; %[Lx1,Ly1 ; Lx2, Ly2; ... ; Lxn, Lyn]
             
             %scatter(obj.Loudspeaker_Locations(:,1),obj.Loudspeaker_Locations(:,2),'r');
             
             width = obj.res * obj.Radius * 2;
             O = width / 2; 
             xyvec = (1:width) - O;
             [xx, yy] = meshgrid(xyvec,xyvec);
             position_vec = [xx(:)'; yy(:)'] / obj.res;
             
             %hold on;
             %scatter(position_vec(1:2:end),position_vec(2:2:end),'g');
             
             p = zeros(width);
             mask = zeros(width);
             rr=obj.Multizone_Soundfield.Radius * obj.res;
             p( (O-rr+1):(O+rr), (O-rr+1):(O+rr) ) = obj.Multizone_Soundfield.Soundfield_desired;
             mask( (O-rr+1):(O+rr), (O-rr+1):(O+rr) ) = obj.Multizone_Soundfield.Soundfield_desired_mask;
             mask=mask(:);
             p = p(:);
             v = p;
             %p = real(p);
             
             %hold off; figure;surf(real(reshape(p,[60 60])),'LineStyle','none');view(2);
             
             %position_vec = [x1, y1;
             %                x2, y2;
             %                xn, yn;]
             position_vec = position_vec(:)';
             %position_vec = [x1, y1, x2, y2, ... , xn, yn;]
             position_vec = repmat(position_vec,[Q0 1]);
             %position_vec = [x1, y1, x2, y2, ... , xn, yn;    ||
             %                x1, y1, x2, y2, ... , xn, yn;    || Q0
             %                x1, y1, x2, y2, ... , xn, yn;]   \/
             position_vec = [reshape(position_vec(:,1:2:end),[1 numel(position_vec)/2])' ...
                             reshape(position_vec(:,2:2:end),[1 numel(position_vec)/2])']; 
             
             %obj.Loudspeaker_Locations = [Lx1, Ly1;
             %                             Lx2, Ly2;
             %                             Lxn, Lyn;]
             spkrPos_vec = repmat(obj.Loudspeaker_Locations,[1 width*width]);
             %spkrPos_vec = [Lx1, Ly1, Lx1, Ly1, ... , Lx1, Ly1;
             %               Lx2, Ly2, Lx2, Ly2, ... , Lx2, Ly2;
             %               Lxn, Lyn, Lxn, Lyn, ... , Lxn, Lyn;]
             spkrPos_vec = [reshape(spkrPos_vec(:,1:2:end),[1 numel(spkrPos_vec)/2])' ...
                            reshape(spkrPos_vec(:,2:2:end),[1 numel(spkrPos_vec)/2])']; 
             
             obj = obj.calc_transferFunc_Pressure(position_vec, spkrPos_vec);
             obj = obj.calc_transferFunc_Velocity(position_vec, spkrPos_vec);
             %obj.transferFunc_Pressure = [Z(xy1,Lxy1);
             %                             Z(xy1,Lxy2);
             %                             ...
             %                             Z(xy1,Lxyn);
             %                             Z(xy2,Lxy1);
             %                             Z(xy2,Lxy2);
             %                             ...
             %                             Z(xy2,Lxyn);
             %                             ...
             %                             Z(xyn,Lxy1);
             %                             Z(xyn,Lxy2);
             %                             ...
             %                             Z(xyn,Lxyn);
             
             obj.transferFunc_Pressure = reshape(obj.transferFunc_Pressure, [Q0 width*width])';
             obj.transferFunc_Velocity = reshape(obj.transferFunc_Velocity, [Q0 width*width])';
             %obj.transferFunc_Pressure = [Z(xy1,Lxy1), Z(xy1,Lxy2), ... , Z(xy1,Lxyn);
             %                             Z(xy2,Lxy1), Z(xy2,Lxy2), ... , Z(xy2,Lxyn);
             %                             ...          ...          ...   ...                             
             %                             Z(xyn,Lxy1), Z(xyn,Lxy2), ... , Z(xyn,Lxyn);
             
             %Pressure
             Zp = obj.transferFunc_Pressure;            
             Zp(mask==0,:)=[];
             p(mask==0,:)=[];
             Beta_p = max(Zp(:)) * 1e-4;
             obj.Loudspeaker_Weights_Pressure = (Zp'*Zp + Beta_p*eye(size(Zp,2))) \ Zp' * p;                       
             
             %Velocity
             Zv = obj.transferFunc_Velocity;            
             Zv(mask==0,:)=[];
             v(mask==0,:)=[];
             Beta_v = max(Zv(:)) * 1e-4;
             obj.Loudspeaker_Weights_Velocity = (Zv'*Zv + Beta_v*eye(size(Zv,2))) \ Zv' * v;                       
        end
        
        function obj = calc_transferFunc_Pressure(obj, P1, P2 )
            
            d = sqrt( sum( (P2 - P1).^2, 2 ) );
            
            omega = 343 * obj.k_global;
            rho_0 = 20*1e-6;
            
            obj.transferFunc_Pressure = 1j * omega * rho_0 * exp(-1j*obj.k_global*d) / 4 / pi ./ d;
        end
        
        function obj = calc_transferFunc_Velocity(obj, P1, P2 )
            
            d = sqrt( sum( (P2 - P1).^2, 2 ) );
            
            k = obj.k_global;
            
            obj.transferFunc_Velocity = 1j*k*exp(-1j*k*d)./(4*pi*d) .* (1 + 1./(1j*k*d));
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
                return;
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
        
        function obj = setRadius(obj, Radius)
           obj.Radius = Radius;           
           obj = obj.calcDesiredMask();
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
            
            obj.Bright_Samples = zeros(O*2,O*2);            
            obj.Bright_Samples_Locations = ones(numel(obj.Bright_Samples), 2)*NaN;
                        
            for x_ = (-O):(O-1)
                for y_ = (-O):(O-1)
                    if obj.Multizone_Soundfield.Bright_Zone.Soundfield_d_mask(x_+O+1,y_+O+1)
                        y = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Bright_Zone.Origin_q.X * obj.res + y_;
                        x = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Bright_Zone.Origin_q.Y * obj.res + x_;
                        
                        Oz = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
                        th2 = -obj.Loudspeaker_Locations(:,1) - pi/2;
                        r2 = obj.Loudspeaker_Locations(:,2) * obj.res;
                        Lp = obj.Loudspeaker_Weights_Pressure;
                        Lv = obj.Loudspeaker_Weights_Velocity;
                        
                        th1 = atan2(Oz-y, Oz-x);
                        r1  = ((Oz-x)^2 + (Oz-y)^2)^(1/2);
                        r = (r1^2 + r2.^2 - 2*r1*r2.*cos(th1-th2)).^(1/2) / obj.res;  %Change to spatial radius (/metre)
                        
                        obj.Bright_Samples_Locations((x_+O)*O*2 + (y_+O+1), :) = [th1, r1 ./ obj.res];
                        
                        omega = 343 * obj.k_global;
                        rho_0 = 20*1e-6;                        
                        Z = 1j * omega * rho_0 * exp(-1j*obj.k_global*r) / 4 / pi ./ r;
                        %H = besselh(0, obj.k_global * r );
                        obj.Bright_Samples(x_+O+1,y_+O+1) = Lp' * Z; % Sum pressures
                    end
                end
            end
            
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
            
            obj.Quiet_Samples = zeros(O*2,O*2);
            obj.Quiet_Samples_Locations = ones(numel(obj.Quiet_Samples), 2) * NaN;
            
            for x_ = (-O):(O-1)
                for y_ = (-O):(O-1)
                    if obj.Multizone_Soundfield.Quiet_Zone.Soundfield_d_mask(x_+O+1,y_+O+1)
                        y = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Quiet_Zone.Origin_q.X * obj.res + y_;
                        x = (obj.res * obj.Radius) + obj.Multizone_Soundfield.Quiet_Zone.Origin_q.Y * obj.res + x_;
                        
                        Oz = length(obj.Soundfield_reproduced) / 2; %Set the centre of the zone for indexing
                        th2 = -obj.Loudspeaker_Locations(:,1) - pi/2;
                        r2 = obj.Loudspeaker_Locations(:,2) * obj.res;
                        Lp = obj.Loudspeaker_Weights_Pressure;
                        Lv = obj.Loudspeaker_Weights_Velocity;
                        
                        th1 = atan2(Oz-y, Oz-x);
                        r1  = ((Oz-x)^2 + (Oz-y)^2)^(1/2);
                        r = (r1^2 + r2.^2 - 2*r1*r2.*cos(th1-th2)).^(1/2) / obj.res;  %Change to spatial radius (/metre)
                                                
                        obj.Quiet_Samples_Locations((x_+O)*O*2 + (y_+O+1), :) = [th1, r1 ./ obj.res];
                        
                        omega = 343 * obj.k_global;
                        rho_0 = 20*1e-6;                        
                        Z = 1j * omega * rho_0 * exp(-1j*obj.k_global*r) / 4 / pi ./ r;
                        %H = besselh(0, obj.k_global * r );
                        obj.Quiet_Samples(x_+O+1,y_+O+1) = Lp' * Z; % Sum pressures
                    end
                end
            end
            
            obj.Quiet_Sample = mean( abs( obj.Quiet_Samples(:) ) );
        end
        
        %% Plotting Functions
        function h = plotSoundfield(obj, field, colour)
            if nargin < 3; colour = jet; end;
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            ZaxisSize = max(abs(real(field(:))));
            f = real(field) .* obj.Desired_Mask;
            CaxisSize = max(abs(f(:))) * 3/4;
            
            XYTick = [1 length(field)/4 length(field)/2 length(field)*3/4 length(field)];
            XYTickLabel = { -length(field)/2/obj.res, -length(field)/4/obj.res, 0, length(field)/4/obj.res, length(field)/2/obj.res};
            
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
            
            obj.zonePlots(real(field));
            obj.speakerPlots(real(field));
            
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
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            % For a Global Zone outline plot
            GlobalZone_temp = obj.Multizone_Soundfield.Bright_Zone;
            GlobalZone_temp.Origin_q = struct('X',0,'Y',0);
            GlobalZone_temp.Radius_q = obj.Multizone_Soundfield.Radius;
            
            maxZ = max(real(field(:)));
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
        
        function speakerPlots(obj, field)
            if nargin < 2; field = obj.Soundfield_reproduced; end;
            
            maxZ = max(real(field(:)));
            O = length(obj.Soundfield_reproduced) / 2 ; %Set the centre of the zone for indexing
            hold on;            
            
            sl_x = obj.Loudspeaker_Locations(:,1) * obj.res;
            sl_y = obj.Loudspeaker_Locations(:,2) * obj.res;

            scatter3(sl_x + O, sl_y + O, ones(length(sl_x),1)*maxZ, 100, 'MarkerEdgeColor', [0 0 0],'MarkerFaceColor', [1 1 1], 'LineWidth', 3);
            
            hold off;
        end
        
    end
    
end

