% The components of the stress tensor MPa
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


%stress tensor
S = [s_xx, s_xy, s_xz; s_xy, s_yy, s_yz; s_xz, s_yz, s_zz];
%normal vector
n = [n_x; n_y; n_z];
%traction vector with respect to xyz axis
t_c = S * n;

%traction with respect to normal and shear
t_n = dot(t_c,n);
t_s = abs(cross(t_c, n));

%princaple stress vectors, and magnitudes
[Evec, Evalue] = eig(S);

disp('traction vector with respect to xyz:');
disp(t_c);
disp('traction normal compontent:');
disp(t_n(1));
disp("shear magnitude: ")
disp(abs(t_s(2)));
disp('principle stress magnitudes: ')
disp("sigma 1: ");
disp(Evalue(1,1));
disp("sigma 2: ");
disp(Evalue(2,2));
disp("sigma 3: ");
disp(Evalue(3,3));

disp("principal stress directions: ")
disp("sigma 1: " );
disp(sEvec(:,1));
disp("sigma 2: ");
disp(sEvec(:,2));
disp("sigma 3: ");
disp(sEvec(:,3));
