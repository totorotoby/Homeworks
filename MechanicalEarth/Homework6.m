

% valley radius (m)
R = 100;
% crustal density (kg/m^3)
rho = 2500;
%gravity (m/s^2)
g = 9.8;

% x coordinates (m)
x = linspace(-400, 400, 1000);
% y coordinates (m)
y = linspace(-300, 0, 1000);
[X, Y] = meshgrid(x,y);
[t, r] = cart2pol(X,Y);
display(t);
display(r);
%stress on radial in radial
display(t);
%display(size(r));
%display(rho*g*R^2 .* 1./r);
Srr = (-1./r)*rho*g*R^2 .* sin(t) + rho*g*r .* sin(t);
%stress on angular in angular
Stt = rho * g * r.* sin(t);

% cut out the valley
Srr(find(r<R)) = nan;
Stt(find(r<R)) = nan;

%contour maps

figure, contourf(X, Y, Srr),  colormap(jet), colorbar
title("sigma rr");

figure, contourf(X, Y, Stt), colormap(jet), colorbar
title("sigma tt");