classdef wideband_soundfield
    % WIDEBAND_SOUNDFIELD - This class represents a wideband (multiple frequency) global
    %                       soundfield containing multiple zones for reproduction by any
    %                       single zone reproduction technique.
    %
    %   Author: Jacob Donley, University of Wollongong, Australia
    %   Email: Jacob.Donley089@uowmail.edu.au
    %
    
    properties        
        res = 20;  % / samples per metre %Resolution of soundfield        
        Soundfield_Wideband = [];
        Error_Field = [];
        Error_Total_MSEdb = 0;
        Soundfield_Desired = [];
        Soundfield_ZoneMask = [];
        Global_Soundfields = [];    % Each global soundfield is for a particular frequency
    end

    methods (Access = private)
            
        function obj = createEmptySoundfield(obj, radius)
            width = int16(obj.res * radius * 2);
            obj.Soundfield_Wideband = zeros(width,width);
        end
        
    end
    
    
    methods
        
        function obj = addGlobalSoundfield(obj, GlobalSoundfield)
            obj.Global_Soundfields = [obj.Global_Soundfields GlobalSoundfield];
        end        
        
        function obj = deleteGlobalSoundfield(obj, index)
            obj.Global_Soundfields = [obj.Global_Soundfields(1:index-1) obj.Global_Soundfields(index+1:end)];
        end
        
        function obj = setResolution(obj, resolution)
           obj.res = resolution;
           for p = 1:length(obj.Global_Soundfields) % For each global soundfield
              obj.Global_Soundfields(p).res = resolution;
              for q = 1:length(obj.Global_Soundfields(p).SpatialZones) % For each spatial zone in that global soundfield
                 obj.Global_Soundfields(p).SpatialZones(q).res = resolution;
              end
           end
        end
        
        function obj = createWidebandSoundfield(obj)   
            
            %Find maximum size radius
            maxRp = 0;
            for p = 1:length(obj.Global_Soundfields)
                if obj.Global_Soundfields(p).getRadius_p_FromZones() > maxRp
                    maxRp = obj.Global_Soundfields(p).getRadius_p_FromZones();
                end
            end
            obj = obj.createEmptySoundfield(maxRp);
            
            % create each global soundfield based on the maximum radius and
            % sum the global soundfields to the wideband field            
            for p = 1:length(obj.Global_Soundfields)
                obj.Global_Soundfields(p) = obj.Global_Soundfields(p).createGlobalSoundfield(maxRp, obj.Soundfield_Wideband);
                obj.Soundfield_Wideband = obj.Global_Soundfields(p).GlobalSoundfield_d;
            end
            
        end
        
        function obj = calcErrorField(obj)
           obj.Error_Field = zeros(length(obj.Soundfield_Wideband), length(obj.Soundfield_Wideband));
           obj.Soundfield_ZoneMask = obj.Error_Field;
           obj.Soundfield_Desired = obj.Soundfield_ZoneMask;
           %obj.Error_Field = obj.Soundfield_Wideband; % Start with the final field
           
           % Subtract each desired spatial zone field from the final field.
           for p = 1:length(obj.Global_Soundfields) % For each global soundfield
              for q = 1:length(obj.Global_Soundfields(p).SpatialZones) % For each spatial zone in that global soundfield
        
                  % Regenerate desired soundfield zone
                  obj.Global_Soundfields(p).SpatialZones(q) = ...
                       obj.Global_Soundfields(p).SpatialZones(q).setDesiredSoundfield('suppress_output');                  
                  desired = obj.Global_Soundfields(p).SpatialZones(q).Soundfield_d .* ...
                      obj.Global_Soundfields(p).SpatialZones(q).Soundfield_d_mask;
                  
                  r = obj.Global_Soundfields(p).SpatialZones(q).Radius_q;
                  xx = obj.Global_Soundfields(p).Radius_p + obj.Global_Soundfields(p).SpatialZones(q).Origin_q.X;
                  xx = round(((xx-r)*obj.res)+1):round(((xx+r)*obj.res));   % TODO: we shoudn't need to round here !
                  yy = obj.Global_Soundfields(p).Radius_p + obj.Global_Soundfields(p).SpatialZones(q).Origin_q.Y;
                  yy = round(((yy-r)*obj.res)+1):round(((yy+r)*obj.res));    % TODO: we shoudn't need to round here !
                                   
                  obj.Soundfield_ZoneMask(yy,xx) = obj.Soundfield_ZoneMask(yy,xx) | ...
                      obj.Global_Soundfields(p).SpatialZones(q).Soundfield_d_mask;
                  
                  obj.Soundfield_Desired(yy,xx) = obj.Soundfield_Desired(yy,xx) + desired;
              end
           end
           
           % Calculate Squared Error field, Mean Squared Error (MSE), Maximum
           % Amplitude from desired field and the MSE in dB with reference
           % to the desired field.
           obj.Error_Field = abs(obj.Soundfield_Desired - obj.Soundfield_Wideband .* obj.Soundfield_ZoneMask) .^ 2;
           
           MSE = mean2(obj.Error_Field);
           MAX = max(abs(obj.Soundfield_Desired(:)));
           obj.Error_Total_MSEdb = 10 * log10(  MSE / MAX^2 );
           
           fprintf('\tError: %fdB \n\n', obj.Error_Total_MSEdb);
           
        end
        
        function plotSoundfield(obj, field, ZAxis_size)
            if nargin < 2; field = obj.Soundfield_Wideband;
                           ZAxis_size = max(abs(real(field(:))));
            end;
            if nargin < 3; ZAxis_size = max(abs(real(field(:))));  
            end; 
            
            subplot(1,2,1);
            surf(real(field),'EdgeColor','None');
            view(2);
            title('Real');
            axis('square');
            axis([1 length(field) 1 length(field)]);
            caxis([-ZAxis_size ZAxis_size]);
            colormap jet; 
            obj.zonePlots();
            
            subplot(1,2,2);
            surf(imag(field),'EdgeColor','None');
            view(2);
            title('Imag');
            axis('square');
            axis([1 length(field) 1 length(field)]);
            caxis([-ZAxis_size ZAxis_size]);
            colormap jet;
            obj.zonePlots();
            
        end
        
        function zonePlots(obj)
            hold on
            th = 0:pi/obj.res:2*pi;
            for p = 1:length(obj.Global_Soundfields)
                for q = 1:length(obj.Global_Soundfields(p).SpatialZones)
                    Ox = obj.Global_Soundfields(p).SpatialZones(q).Origin_q.X * obj.res + length(obj.Global_Soundfields(p).GlobalSoundfield_d) / 2 + 0.5;
                    Oy = obj.Global_Soundfields(p).SpatialZones(q).Origin_q.Y * obj.res + length(obj.Global_Soundfields(p).GlobalSoundfield_d) / 2 + 0.5;
                    xunit = obj.Global_Soundfields(p).SpatialZones(q).Radius_q * obj.res * cos(th) + Ox;
                    yunit = obj.Global_Soundfields(p).SpatialZones(q).Radius_q * obj.res * sin(th) + Oy;
                    plot(xunit, yunit);
                end
            end
            hold off
        end
        
        
        function plotError(obj)
            ZaxisSize = max(abs(real(obj.Soundfield_Desired(:))));
            
            %subplot(1,2,1);
            surf(real(obj.Error_Field),'EdgeColor','None');
            title('Real');
            axis('square');
            axis([1 length(obj.Error_Field) 1 length(obj.Error_Field) -ZaxisSize ZaxisSize]);
            view([0 0]);
            caxis([-ZaxisSize ZaxisSize]);
            colormap jet; 
            obj.zonePlots();
            
%             subplot(1,2,2);
%             surf(imag(obj.Error_Field),'EdgeColor','None');
%             title('Imag');
%             axis('square');
%             axis([1 length(obj.Error_Field) 1 length(obj.Error_Field) -ZaxisSize ZaxisSize]);            
%             view([0 0]);
%             %caxis([-ZaxisSize ZaxisSize]);
%             colormap jet;
%             obj.zonePlots();
%             %obj.plotSoundfield(obj.Soundfield_Desired);
        end
    end
    
end

