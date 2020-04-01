% elliptical_hole.m
% Calculate stress and displacement fields in an infinite
% elastic body with an elliptcal hole following:
% Pollard, 1973, Mathematical Geology, Vol.5, No.1, p. 11-25.
% Boundary conditions are a biaxial remote stress and internal pressure. 
% The elastic body is linear, homogeneous, and isotropic.
% Stress and displacement components are positive as depicted in Figure 2.
% That is compression is positive and displacement components are positive
% in the negative coordinate directions.
% See also Mushkelishvili 1953, and Stevenson, 1945
% M-file by JOK 21.3.2005 modified by DDP 26.07.2005

clear all % clear memory and figures
close all
% input geometry, boundary conditions, elastic moduli
a = 1; % semi-major axis length 
b = .9; % semi-minor axis length, b < a
beta = 0; % angle (degrees) of s3 to x-axis
s3 = -1; % remote stress acting at angle beta to x-axis
s1 = -1; % remote stress acting at angle beta + pi/2 to x-axis
p = 2; % pressure in the hole
pr = 0.25; % Poisson's ratio
g = 1e3; % elastic shear modulus

x=-2*a:.005*a:2*a;   % define x-coords. of grid
y=-2*a:.005*a:2*a;   y=y+1e-6; % define y-coords. of grid
[X,Y] = meshgrid(x,y); % define Cartesian grid

% elliptical coordinates for the grid
f = sqrt(a^2-b^2); % focal length of elliptical hole
XI = acosh((sqrt((X+f).^2+Y.^2)+sqrt((X-f).^2+Y.^2))/(2*f));
ET = atan2(Y,(X.*tanh(XI)));

xio = atanh(b/a); % xi coordinate of hole boundary
ET(XI<xio) = nan; % set coords within hole to nan
XI(XI<xio) = nan; % set coords within hole to nan

% plot elliptical coordinate curves and hole boundary
% contour(X,Y,XI,25,'k'), hold on 
% contour(X,Y,ET,25,'k'), axis image
% XO = f*cosh(xio)*cos(0:pi/50:2*pi); 
% YO = f*sinh(xio)*sin(0:pi/50:2*pi);
% plot(XO,YO,'b-'), xlabel('x/|a|'), ylabel('y/|a|')
% title('elliptical coordinate system and hole')
% 

p = -p;

% shorthands
S = [s3,s1]; % vector of remote principal stresses
s2xio = sinh(2*xio); c2xio = cosh(2*xio);
e2xio = s2xio + c2xio;

S2XI = sinh(2*XI); S2ET = sin(2*ET);
DEN = cosh(2*XI)-cos(2*ET); DEN3 = DEN.^3;

I = sinh(XI).*sin(ET);
J = cosh(XI).*cos(ET);
K = sinh(XI).*cos(ET);
L = cosh(XI).*sin(ET);

A = J.*K.*K.*K-3*J.*K.*L.*L-3*I.*K.*K.*L+I.*L.*L.*L;
B = I.*K.*K.*K-3*I.*K.*L.*L+3*J.*K.*K.*L-J.*L.*L.*L; MB = -B;
C = J.*K.*K.*K-3*J.*K.*L.*L+3*I.*K.*K.*L-I.*L.*L.*L;
D = -I.*K.*K.*K+3*I.*K.*L.*L+3*J.*K.*K.*L-J.*L.*L.*L;
E = I.*K.*K.*K-3*J.*K.*K.*L-3*I.*K.*L.*L+J.*L.*L.*L;

% Stress and displacement components for 
% solution (3) uniform hydrostatic pressure p everywhere
SXX = p*ones(size(X)); SYY = SXX; SXY = 0*SXX;
UX = .25*p*(f/g)*(J*(3-4*pr)-J);
UY = .25*p*(f/g)*(I*(3-4*pr)-I);

% Loop to compute stress and displacement components for
% solution (1) remote uniaxial stress s3-p at angle beta to x-axis 
% solution (2) remote uniaxial stress s1-p at angle beta+pi/2
for k = 1:2
    be = beta*pi/180 + (k-1)*pi/2;
    s2b = sin(2*be); c2b = cos(2*be);
    smp = S(k)-p;
      
% Derivatives of the stress functions
REPSIP = .25*smp*(e2xio*c2b+(S2XI-e2xio*(S2XI*c2b+S2ET*s2b))./DEN);

REZPSIPP = -2*smp*(((1-e2xio*c2b))*A-e2xio*s2b*B)./DEN3;         

IMZPSIPP = -2*smp*(((1-e2xio*c2b))*MB-e2xio*s2b*A)./DEN3;

RECHIPP = -2*smp*(.25*e2xio*(c2xio*c2b-(s2xio*c2b*S2XI+c2xio*s2b*S2ET)./DEN)...
    +((c2b-c2xio+e2xio*s2xio*c2b)*C+e2xio*c2xio*s2b*D)./DEN3);

IMCHIPP = -2*smp*(.25*e2xio*(s2xio*s2b-(c2xio*s2b*S2XI-s2xio*c2b*S2ET)./DEN)...
    +((c2b-c2xio+e2xio*s2xio*c2b)*E+e2xio*c2xio*s2b*C)./DEN3);
                 
REPSI = .25*smp*f*(e2xio*((J-K)*c2b+L*s2b)+K);

IMPSI = .25*smp*f*(e2xio*((I-L)*c2b-K*s2b)+L);

REZPSIP = .25*smp*f*(e2xio*J*c2b+(J.*S2XI-I.*S2ET)./DEN...
    -e2xio*((S2XI.*(J*c2b+I*s2b)-S2ET.*(I*c2b-J*s2b))./DEN));

IMZPSIP = .25*smp*f*(e2xio*I*c2b+(I.*S2XI-J.*S2ET)./DEN...
    -e2xio*((S2XI.*(I*c2b+J*s2b)-S2ET.*(J*c2b-I*s2b))./DEN));

RECHIP = -.5*smp*f*((K*(c2xio-c2b)-e2xio*(K*s2xio*c2b+L*c2xio*s2b))./DEN...
    +e2xio*(c2b*(J*c2xio-K*s2xio)-s2b*(I*s2xio-L*c2xio)));

IMCHIP = -.5*smp*f*((L*(c2xio-c2b)-e2xio*(L*s2xio*c2b+K*c2xio*s2b))./DEN...
    +e2xio*(c2b*(-I*c2xio+L*s2xio)-s2b*(J*s2xio-K*c2xio)));

% Stress components
SXX = SXX + 2*REPSIP - REZPSIPP - RECHIPP;
SYY = SYY + 2*REPSIP + REZPSIPP + RECHIPP;
SXY = SXY + IMZPSIPP + IMCHIPP;

% Displacement components
UX = UX + (1/2*g)*((3-4*pr)*REPSI - REZPSIP - RECHIP);
UY = UY + (1/2*g)*((3-4*pr)*IMPSI - IMZPSIP - IMCHIP);

end

% % % Eliminate points with extreme values of stress
% SXXm = max(max(SXX));
% SYYm = max(max(SYY));
% SXYm = max(max(SXY));
% X(SXX > (0.75*SXXm)) = nan;  
% X(SYY > (0.75*SYYm)) = nan;  
% X(SXY > (0.75*SXYm)) = nan;  

% eliminate the singularity
SYY(abs(X)<=a&abs(Y)<=.01) = nan;
SXY(abs(X)<=a&abs(Y)<=.01) = nan;
SXX(abs(X)<=a&abs(Y)<=.01) = nan;

% Plot stress components
figure, contourf(X,Y,(SXX),15), axis image,  colorbar
title('\sigma_{xx}'), ylabel('y/|a|'), xlabel('x/|a|')
figure, contourf(X,Y,(SYY),15), axis image, colorbar
title('\sigma_{yy}'), ylabel('y/|a|'), xlabel('x/|a|')
figure, contourf(X,Y,(SXY),15), axis image, colorbar
title('\sigma_{xy}'), ylabel('y/|a|'), xlabel('x/|a|')


    % Calculate principal stresses
    magS1 = [];
    magS2 = [];
    for i = 1:length(x)
        for j = 1:length(x)
            
            if isnan(SXX(i,j))
                magS1(i,j) = nan;
                magS2(i,j) = nan;
                
            else
            
            ST = [SYY(i, j) SXY(i, j); SXY(i,j) SXX(i,j)];
            [S,V] = eig(ST);          
                
                magS1(i,j) = V(1,1);
                magS2(i,j) = V(2,2);
                
                
                
            end
           
        end
    end
    
    % eliminate the singularity
magS1(abs(X)<=a&abs(Y)<=.001) = nan;
magS2(abs(X)<=a&abs(Y)<=.001) = nan;

    
    % close all
    % % plot stress components as contour maps
%     figure, contourf(X/a,Y/a,magS1/(p-s1),20), axis equal tight
%     ylabel('y/a','fontsize',18,'fontweight','bold','fontname','Times New Roman');
%     xlabel('x/a','fontsize',18,'fontweight','bold','fontname','Times New Roman');
%     title('\sigma_{1}/\Delta \sigma_{I}', 'fontsize',18,'fontweight','bold','fontname','Times New Roman'), colorbar
%     set(gca,'fontsize',18)

%     figure, contourf(X,Y,magS1,20),  axis image,  colorbar
%     ylabel('y/a','fontsize',18,'fontweight','bold','fontname','Times New Roman');
%     xlabel('x/a','fontsize',18,'fontweight','bold','fontname','Times New Roman');
%     title('\sigma_{1}', 'fontsize',18,'fontweight','bold','fontname','Times New Roman'), colorbar
%     set(gca,'fontsize',18)
% 
