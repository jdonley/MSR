function [ k ] = getAliasingFrequency( setup, method )
%GETALIASINGFREQUENCY Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    method = 'new';
end

%% Known Constants
width = double(int16(setup.res * setup.Radius * 2));
O = [width / 2; ...                                             %Origin
    width / 2; ...
    0 * setup.res];
B =[setup.Multizone_Soundfield.Bright_Zone.Origin_q.X;...       % Bright zone center point
    setup.Multizone_Soundfield.Bright_Zone.Origin_q.Y; ...
    0] * setup.res + O;
[Bth,Br]=cart2pol(B(1)/setup.res, B(2)/setup.res);              % Bright zone angle and distance
[Bth_,Br_]=cart2pol((B(1)-O(1))/setup.res, (B(2)-O(2))/setup.res);    % Bright zone angle and distance from origin
Q =[setup.Multizone_Soundfield.Quiet_Zone.Origin_q.X;...        % Quiet zone center point
    setup.Multizone_Soundfield.Quiet_Zone.Origin_q.Y; ...
    0] * setup.res + O;
[Qth,Qr]=cart2pol(Q(1)/setup.res, Q(2)/setup.res);              % Quiet zone angle and distance
[Qth_,Qr_]=cart2pol((Q(1)-O(1))/setup.res, (Q(2)-O(2))/setup.res);    % Quiet zone angle and distance from origin
r1 = setup.Multizone_Soundfield.Bright_Zone.Radius_q;           % Bright zone radius
r2 = setup.Multizone_Soundfield.Quiet_Zone.Radius_q;            % Quiet zone radius
theta = setup.Multizone_Soundfield.Bright_Zone.SourceOrigin.Angle/180*pi;   % Angle of bright zone desired wave
L = setup.Loudspeaker_Count;                                    % Number of loudspeakers
D_l = setup.Speaker_Array_Length;                                 % Length of loudspeaker array (for linear array)
phi_c = setup.Speaker_Array_Centre/180*pi;                        % Angle to center of array
Rl = (setup.Radius);                                            % Loudspeaker radius (for circular array) OR Perpendicular distance to center of linear array
phi_L=setup.Speaker_Arc_Angle/180*pi;                           % Angle of loudspeaker array (angle from begining to end)
R = setup.Multizone_Soundfield.Radius;
Rprime = max([r1+Br_,r2+Qr_]);                                    % Smallest circular radius encompassing all zones

% Some common vectors
Bv = B-O; % Vector for point B
Qv = Q-O; % Vector for point Q

%% Compute Unknown Values
if strcmpi(setup.Speaker_Array_Type(1),'2')
    setup.Speaker_Array_Type(1) = [];
    L = L/2;
end
if strcmpi('line',setup.Speaker_Array_Type) ...
        || strcmpi('plane',setup.Speaker_Array_Type)
    %% Planar Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take line from planar array and assume aliasing frequency is the same
    if strcmpi('plane',setup.Speaker_Array_Type)
       L = sqrt(L); 
    end
    
    %% Linear Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = O + setup.res*Rl*[...                               % Center point of loudspeaker array
        cos(phi_c) -sin(phi_c) 0;...      
        sin(phi_c)  cos(phi_c) 0;...
        0           0          1]*[1;0;0];
    m1 = [cos(theta);sin(theta);0];                         % Planewave unit vector
    m2 = [cos(phi_c-pi/2);sin(phi_c-pi/2);0];               % Loudspeaker array unit vector
    A_=[m1 -m2];	b_ = [C-B]; w = pinv(A_)*b_;            % Solve intersection of planewave and loudspeaker array
    P = B+w(1)*m1; P_ = C+w(2)*m2;                          % Find point of intersection
    if ~all(round(P,10)==round(P_,10))                      % Check if intersection occurred
        error('Planewave does not intersect loudspeaker array.');
    end
    
    Pv = P-O;                                               % Intersection point vector
    PQv = Qv-Pv;                                            % Intersection to Quiet zone vector
    d = sqrt(sum((Q-P).^2))/setup.res;                      % Length of PQ vector
    gam_found=[];
    for pm = [-1 1]                                         % Find tangent either side of the zone
        if abs(r1+r2) <= d
        alph_=pm*asin((r1+r2)/d);                           % Maximum allowable angle of grating lobe
                                                            % Numerator in arcsin gives distance from quiet zone center to consider aliasing frequency
        else alph_=pm*pi/2;
        end
        PAv=[cos(alph_) -sin(alph_) 0;...                   % Tangent
             sin(alph_)  cos(alph_) 0;...
             0           0          1] * PQv;                        
        PBv = Bv-Pv;                                        % Intersection to Bright zone vector
        gam_found(end+1) = atan2(norm(cross(PAv,PBv)), dot(PAv,PBv));   % One of two maximum allowable angles for the grating lobe
        % MATHEMATICALLY EQUIVALENT --> acos(dot(v1, v2) / (norm(v1) * norm(v2)))        
    end
    
    gam_found = max(gam_found);                             % Choose the angle farthest from the planewave vector
    Theta = abs(pi - abs(theta) - abs(phi_c) );             % Find adjusted angle from centre of array
    k =   2*pi*(L-1) / (abs(sin(gam_found - Theta) + sin(Theta))*D_l) ;    % Equivalent grating lobe frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
elseif strcmpi('circle',setup.Speaker_Array_Type)
    %% Circular Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perpendicular distance from planewave vector to center
    d_pb = abs( Br_*sin(-theta + Bth_) );
    
    % Angle to point of intersection of planewave vector and loudspeaker
    % array
    alpha_ = pi +theta - asin(d_pb/Rl);
    
    % Find tangent
    P = O + Rl*setup.res * ...
        [cos(alpha_) -sin(alpha_) 0;...
        sin(alpha_)  cos(alpha_) 0;...
        0           0       1]*[1;0;0];    
    Pv = P-O;
    PQv = Qv-Pv;
    d = sqrt(sum((Q-P).^2))/setup.res;
    
    d_u=[];
    for pm = [-1 1]
        %alph_ is angle of one of two tangents
        alph_ = pm * asin((r1+r2)/d);
        
        %PAv is one of two tangents
        PAv=[cos(alph_) -sin(alph_) 0;...
            sin(alph_)  cos(alph_) 0;...
            0           0      1] * PQv;
        
        rotMat_90 = [cos(pi/2) -sin(pi/2) 0;...
            sin(pi/2)  cos(pi/2) 0;...
            0           0      1];
        % d_u is perpendicular distance of one tangent to the center
        d_u(end+1) = abs(  Pv.'   * (rotMat_90*PAv  )) ...
            / sqrt(sum((PAv  ).^2)) / setup.res;
    end
    d_u = max(d_u); %Choose the farthest from the center (outermost tangent)
    
    
    % Maximum allowable aliasing frequency
    if strcmpi('pw',setup.Multizone_Soundfield.Bright_Zone.SourceType) && strcmpi(method,'new')
        k = max( ...
            [(2*pi*L - phi_L) / ((d_u + d_pb)*phi_L), ...
            (2*pi*L - phi_L) / (2*Rprime*phi_L)]);
    elseif strcmpi(method,'old')
        k = (2*pi*L - phi_L) / (2*Rprime*phi_L);
    elseif strcmpi(method,'old2')
        k = (2*pi*L - phi_L) / (2*R*phi_L);
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

end

