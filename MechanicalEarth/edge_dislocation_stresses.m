% plot Cartesian stress components near the edge dislocation
% Weertman & Weertman (1964)

%clear all, clf reset; % clear memory and figures
b = .01; % b is Burgers vector
mu = 30000; pr = 0.25; lm = (2*mu*pr)/(1-2*pr); % Elastic moduli
c1 = b/(2*pi); c2 = 2*mu*(mu+lm)/(2*mu+lm);
x = linspace(-1,1,100)+eps;
y = linspace(-1,1,100);
[X,Y]=meshgrid(x,y); % define Cartesian grid
DEN = (X.^2 + Y.^2).^2;
SXX = c1*c2*((Y.*(3*X.^2 + Y.^2))./DEN);
SYY = -c1*c2*((Y.*(X.^2 - Y.^2))./DEN);
SXY = -c1*c2*((X.*(X.^2 - Y.^2))./DEN);
[T,R] = cart2pol(X,Y); % R is radial distance from dislocation
%SXX(find(R<(5*b))) = nan; % avoid dislocation core (R < 5b)
SYY(find(R<(5*b))) = nan; SXY(find(R<(5*b))) = nan;

% plot contour maps
figure, contourf(X./b,Y./b,SXX,25), axis equal ij tight, colormap(jet)
title('stress \sigma_{xx}'), xlabel('x/b'), ylabel('y/b'), colorbar
figure, contourf(X./b,Y./b,SYY,25), axis equal ij tight, colormap(jet)
title('stress \sigma_{yy}'), xlabel('x/b'), ylabel('y/b'), colorbar
figure, contourf(X./b,Y./b,SXY,25), axis equal ij tight, colormap(jet)
title('stress \sigma_{xy}'), xlabel('x/b'), ylabel('y/b'), colorbar