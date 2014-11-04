%function [H, G, L, g] = calculate_jacobians()

syms x y z phi theta psi Vx Vy Vz Wx Wy Wz wx wy wz dt ex ey ez ephi etheta epsi

h  = [x y z phi theta psi]';
xt = [x y z phi theta psi];
H  = jacobian(h,xt);
g1 = matlabFunction(H, 'file', 'H_function', 'vars', {[x; y; z], [phi; theta; psi]});

% R taken from Project one.
R = [cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta), -cos(phi)*sin(theta); 
                                     -cos(phi)*sin(psi),                                  cos(phi)*cos(psi),             sin(phi);
     cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi), sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi), cos(phi)*cos(theta)];
 
W_SkSym = [           0, -(Wz+epsi),  Wy+etheta;
                Wz+epsi,          0, -(Wx+ephi);
           -(Wy+etheta),    Wx+ephi,          0];
       
%R_dot = R * W_SkSym;
Delta_R = eye(3) + W_SkSym * dt;
[phin,thetan,psin] = RotToRPY_ZXY_sym(R * Delta_R);
to_eul = [phin,thetan,psin]';
noise = [ex ey ez ephi etheta epsi];

%%
g = [[x y z]' + R * ([Vx Vy Vz] + [ex ey ez])' * dt; to_eul];
G = jacobian(g,xt);
g2 = matlabFunction(G, 'file', 'G_function', 'vars', {[x; y; z; phi; theta; psi], [Vx; Vy; Vz; Wx; Wy; Wz], [ex; ey; ez; ephi; etheta; epsi], dt});

%
L = jacobian(g,noise);
g3 = matlabFunction(L, 'file', 'L_function', 'vars', {[x; y; z; phi; theta; psi], [Vx; Vy; Vz; Wx; Wy; Wz], [ex; ey; ez; ephi; etheta; epsi], dt});
g4 = matlabFunction(g, 'file', 'g_func', 'vars', {[x; y; z; phi; theta; psi], [Vx; Vy; Vz; Wx; Wy; Wz], [ex; ey; ez; ephi; etheta; epsi], dt});

%end

%% Old Code

% wx = Wx/norm([Wx Wy Wz]);
% wy = Wy/norm([Wx Wy Wz]);
% wz = Wz/norm([Wx Wy Wz]);
% wt = dt*norm([Wx Wy Wz]);
% 
% c = cos(wt);
% s = sin(wt);
% C = 1 - c;
% 
% dR = [     wx * wx * C + c, wx * wy * C - wz * s, wx * wz * C + wy * s;
%       wy * wx * C + wz * s,      wy * wy * C + c, wy * wz * C - wx * s;
%       wz * wx * C - wy * s, wz * wy * C + wx * s,      wz * wz * C + c];
% 