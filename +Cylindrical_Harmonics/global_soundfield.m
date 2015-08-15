classdef global_soundfield
    % GLOBAL_SOUNDFIELD - This class represents a global soundfield
    %                     containing multiple zones for reproduction by any
    %                     single zone reproduction technique.
    %
    %   This class was implemented based on the publication:
    %   Y. J. Wu and T. D. Abhayapala, "Spatial multizone soundfield
    %   reproduction: Theory and design," Audio, Speech, and Language
    %   Processing, IEEE Transactions on, vol. 19, pp. 1711-1720, 2011.
    %
    %   Author: Jacob Donley, University of Wollongong, Australia
    %   Email: Jacob.Donley089@uowmail.edu.au
    %
    
    properties        
        res = 20;  % / samples per metre %Resolution of soundfield
        Radius_p = 0.0;
        k_global = 0;
        Mode_forced = -1;
        SpatialZones;
        GlobalSoundfield_d = []; %The complex values of the global desired sound field.
        Beta_Coeffs = [];           % Beta_Coeffs are the global soundfield coefficients for all zones
    end
    
    

    methods (Access = private)
            
        function obj = createEmptyGlobalSoundfield_d(obj)
            width = int16(obj.res * obj.Radius_p * 2);
            obj.GlobalSoundfield_d = zeros(width,width);
        end
        
    end
    
    
    methods
        function obj = addSpatialZone(obj, spatialZone, Distance, Angle_degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Check for overlapping soundfields here before adding it %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(obj.SpatialZones) < 1
                obj.k_global = spatialZone.getWavenumber();
            elseif obj.k_global ~= spatialZone.getWavenumber()
                error(['This spatial zone does not match the frequency/wavenumber of the global soundfield it is being added to.' ...
                      '\n Global Soundfield Frequency = %0.2fHz\tGlobal Soundfield Wavenumber = %0.2f' ...
                      '\n Spatial Zone Frequency = %0.2fHz\t\tSpatial Zone Wavenumber = %0.2f'], ...
                      obj.getFrequency(), obj.k_global, spatialZone.Frequency, spatialZone.getWavenumber());
                return;
            end
            Angle_rads = Angle_degrees * 2*pi / 360 ;
            [X, Y] = pol2cart(Angle_rads, Distance);
            spatialZone = spatialZone.setOrigin(X, Y);
            obj.SpatialZones = [obj.SpatialZones spatialZone];
        end        
        function obj = deleteSpatialZone(obj, index)
            obj.SpatialZones = [obj.SpatialZones(1:index-1) obj.SpatialZones(index+1:end)];
        end 
        function f = getFrequency(obj)
            f = obj.k_global * 343 / (2 * pi);
        end
        function obj = calc_Beta_Coeffs(obj)
            Q = length(obj.SpatialZones); % Total number of zones
            M0 = obj.getGlobalModeLimit();
            Alpha = [];
            T = [];
            for q = 1:Q         % Create Q simultaneous equations               
                Alpha = [Alpha obj.SpatialZones(q).Alpha_Coeffs];  % Create Alpha matrix
                
                Tq = [];
                Mq = fix( obj.SpatialZones(q).getWavenumber * exp(1) * obj.SpatialZones(q).Radius_q / 2);
                for m = (-Mq:Mq)
                    T_hrz = [];
                    for n = (m+M0):-1:(m-M0)
                        [r_q0, theta_0q] = obj.SpatialZones(q).getOriginInPolarCoords();
                        T_0q = besselj(n, obj.SpatialZones(q).getWavenumber * r_q0 ) * exp( -1j * n * (theta_0q+pi));
                        T_hrz = [T_hrz T_0q];             % Create T matrix
                    end
                    Tq = [Tq; T_hrz];
                end
                if q==1
                T = [T; Tq*100];
                else
                T = [T; Tq*0.05];                    
                end
            end
            Alpha = transpose(Alpha);
            obj.Beta_Coeffs = transpose(pinv(T) * Alpha);
        end
        
        function M0 = getGlobalModeLimit(obj)
            M0 = fix( obj.k_global * exp(1) * obj.Radius_p / 2); 
            if obj.Mode_forced ~= -1; M0 = obj.Mode_forced; end;
        end
        
        function radius_p = getRadius_p_FromZones(obj)
           Q = length(obj.SpatialZones); % Total number of zones
           
           % Find the radius of the smallest circle that encloses all spatial zones of interest.
           radius_p = obj.SpatialZones(1).Radius_q + obj.SpatialZones(1).getOriginInPolarCoords();
            for q = 2:Q 
                if (obj.SpatialZones(q).Radius_q + obj.SpatialZones(q).getOriginInPolarCoords()) > radius_p
                   radius_p = obj.SpatialZones(q).Radius_q + obj.SpatialZones(q).getOriginInPolarCoords();
                end
            end
        end
        
        function obj = createGlobalSoundfield(obj, Radius, Existing_Field)
            if nargin < 2
                obj.Radius_p = obj.getRadius_p_FromZones;
            else
               obj.Radius_p = Radius; 
            end
            M0 = obj.getGlobalModeLimit; 
            
            obj = obj.calc_Beta_Coeffs();
            if nargin < 3
                obj = obj.createEmptyGlobalSoundfield_d;
            else
                obj.GlobalSoundfield_d = Existing_Field;
            end
            O = length(obj.GlobalSoundfield_d) / 2; %Set the centre of the zone for indexing
                        
            n=0; 
            fprintf('--- Creating Global Soundfield ---\n');
            fprintf('\t%d zones at %.0fHz.\n\tCompletion: ',length(obj.SpatialZones), obj.getFrequency);            
            
            hold on
            for m = -M0:M0
                for x = 1:O*2
                    for y = 1:O*2 
                        [theta, r] = cart2pol(x - O, y - O); %Convert to polar coords from centre of the zone
                        r = r / obj.res;  %Change to spatial radius (/metre)                                    
                        J = besselj(m, obj.k_global * r );
                        e = exp(1j * m * theta);
                        obj.GlobalSoundfield_d(y, x) = obj.GlobalSoundfield_d(y, x) + obj.Beta_Coeffs(m+M0+1) * J * e; % Sum pressures
                    end
                    if round(mod(x/(O*2)*100,20))==0
                        fprintf(repmat('\b',1,n));
                        n=fprintf('%.2f%%', (x + O*2*(m+M0)) / (O*2*(2*M0+1)) * 100);
                    end
                end
                obj.plotSoundfield();
                drawnow  
            end
            hold off
            fprintf('\n\n');
        end
        
        function plotSoundfield(obj)
            ZaxisSize = max(abs(real(obj.GlobalSoundfield_d(:))));
            
            subplot(1,2,1);
            surf(real(obj.GlobalSoundfield_d),'EdgeColor','None');
            view(2);
            title('Real');
            axis('square');
            axis([1 length(obj.GlobalSoundfield_d) 1 length(obj.GlobalSoundfield_d)]);
            caxis([-ZaxisSize ZaxisSize]);
            colormap gray; 
            obj.zonePlots();
            
            subplot(1,2,2);
            surf(imag(obj.GlobalSoundfield_d),'EdgeColor','None');
            view(2);
            title('Imag');
            axis('square');
            axis([1 length(obj.GlobalSoundfield_d) 1 length(obj.GlobalSoundfield_d)]);
            caxis([-ZaxisSize ZaxisSize]);
            colormap jet;
            obj.zonePlots();
            
%             subplot(2,2,3);
%             surf(real(obj.GlobalSoundfield_d));
%             title('Real');
%             axis('square');
%             subplot(2,2,4);
%             surf(imag(obj.GlobalSoundfield_d));
%             title('Imag');
%             axis('square');
        end
        
        function zonePlots(obj)
            hold on
            th = 0:pi/obj.res:2*pi;
            for q = 1:length(obj.SpatialZones)
                Ox = obj.SpatialZones(q).Origin_q.X * obj.res + length(obj.GlobalSoundfield_d) / 2 + 1;
                Oy = obj.SpatialZones(q).Origin_q.Y * obj.res + length(obj.GlobalSoundfield_d) / 2 + 1;
                xunit = obj.SpatialZones(q).Radius_q * obj.res * cos(th) + Ox;
                yunit = obj.SpatialZones(q).Radius_q * obj.res * sin(th) + Oy;
                plot(xunit, yunit);
            end
            hold off
        end
        
    end
    
end

