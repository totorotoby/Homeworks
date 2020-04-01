% The components of the stress tensor
s_xx = -40;
s_yy = -60;
s_zz = -80;
s_xy = 20;
s_xz = -10;
s_yz = 10;
% The components of the normal vector
n_x = .3015;
n_y = .3015;
n_z = .9045;

%angle from the x-axis
alpha = atan2(sqrt(n_z^2+n_y^2), n_x);


%stress tensor
S = [s_xx, s_xy, s_xz; s_xy, s_yy, s_yz; s_xz, s_yz, s_zz];
%normal vector
n = [n_x; n_y; n_z];
%traction vector with respect to xyz axis
t_c = S * n;

%traction with respect to normal and shear
t_ns = [s_xx*cos(alpha)^2 + s_yy*cos(alpha)^2 + 2*s_xy*sin(alpha)*cos(alpha);
    -(s_xx - s_yy) * sin(alpha) * cos(alpha) + s_xy*(cos(alpha)^2 - sin(alpha)^2)];


%princaple stress vectors, and magnitudes
[Evec, Evalue] = eig(S);

disp('traction vector with respect to xyz:');
disp(t_c);
disp('traction normal compontent:');
disp(t_ns(1));
disp("shear magnitude: ")
disp(abs(t_ns(2)));
disp('principle stress magnitudes: ')
disp(Evalue(1,1));
disp(Evalue(2,2));
disp(Evalue(3,3));

disp("principal stress directions: ")
disp(Evec(:,1));
disp(Evec(:,2));
disp(Evec(:,3));
