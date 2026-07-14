function res = f_dyn_1(t,phi,v,b,phi_0,g,f,za,kh,M)
% Definition of 2nd order differential equation
if abs(phi - pi/2) < 1e-6
        phi = pi/2 + 1e-6; % Adjust phi to avoid cos(phi) being zero
    elseif abs(phi + pi/2) < 1e-6
        phi = -pi/2 - 1e-6; % Adjust phi to avoid cos(phi) being zero
    end

    res = 1/(2*b*cos(phi)) * ( 2*b*sin(phi)*v.^2 +4*pi^2*f.^2*za*sin(2*pi.*f.*t) ...
        -g  + (9*kh*b)/(2*M)*tan(phi).*(cos(phi) - cos(phi_0)  )       );
    
%% phi limitation
% if abs(phi - pi/2) < 1e-6 || abs(phi + pi/2) < 1e-6
%         res = 0; % Treat as rigid body at boundaries
%     else
%         res = 1/(2*b*cos(phi)) * ( 2*b*sin(phi)*v.^2 + 4*pi^2*f.^2*za*sin(2*pi.*f.*t) ...
%             -g  + (9*kh*b)/(2*M)*tan(phi).*(cos(phi) - cos(phi_0)) );
%     end
% 


end