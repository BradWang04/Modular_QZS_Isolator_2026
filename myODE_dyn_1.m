function dphidt = myODE_dyn_1(t,phi,b,phi_0,g,f,za,kh,M)

    dphidt = [phi(2);           % phi' = v
              f_dyn_1(t,phi(1),phi(2),b,phi_0,g,f,za,kh,M)];% v' = f(t,phi,v)

%% phi limitation

% Ensure phi stays within the range [-pi/2, pi/2]
%     if phi(1) > pi/2
%         phi(1) = pi/2;
%         phi(2) = 0; % Set velocity to 0 at the boundary
%     elseif phi(1) < -pi/2
%         phi(1) = -pi/2;
%         phi(2) = 0; % Set velocity to 0 at the boundary
%     end
%     dphidt = [phi(2);           % phi' = v
%               f_dyn_1(t,phi(1),phi(2),b,phi_0,g,f,za,kh,M)];% v' = f(t,phi,v)

end