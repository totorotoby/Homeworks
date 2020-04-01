% plot stress components for circular hole in infinite plate
% biaxial remote stress
% Jaeger and Cook (1979)
% equations 6.108 - 6.110 in Pollard and Fletcher (2005)

clear all, clf reset; close all % clear memory and figures
sH = 1; sh = 1; pm = 2;
x = linspace(-4,4,161)+eps; y = linspace(-4,4,161);
[X,Y] = meshgrid(x,y);
[TH,R] = cart2pol(X,Y);
ST = sin(TH); S2T = sin(2*TH); ST2 = ST.^2; 
CT = cos(TH); C2T = cos(2*TH); CT2 = CT.^2;
ri = 1; 
R2 = (ri./R).^2; R4 = R2.^2;

% Polar stress components
SRR = -(0.5*(sH+sh)*(1-R2))-(pm*R2)-(0.5*(sH-sh)*((1-4*R2+3*R4).*C2T));
STT = -(0.5*(sH+sh)*(1+R2))+(pm*R2)+(0.5*(sH-sh)*((1+3*R4).*C2T));
SRT = 0.5*(sH-sh)*((1+2*R2-3*R4).*S2T);
SRR(find(R<1)) = nan; STT(find(R<1)) = nan; SRT(find(R<1)) = nan;
subplot(2,2,1), contourf(X,Y,SRR,10), colorbar, title('stress srr')
subplot(2,2,2), contourf(X,Y,STT,10), colorbar, title('stress stt')
subplot(2,2,3), contourf(X,Y,SRT,10), colorbar, title('stress srt')

% Cartesian stress components
SXX = SRR.*CT2+STT.*ST2-2*SRT.*CT.*ST;
SYY = SRR.*ST2+STT.*CT2+2*SRT.*CT.*ST;
SXY = (SRR-STT).*CT.*ST+SRT.*(CT2-ST2);
figure, subplot(2,2,1), contourf(X,Y,SXX,10), colorbar, title('stress sxx')
subplot(2,2,2), contourf(X,Y,SYY,10), colorbar, title('stress syy')
subplot(2,2,3), contourf(X,Y,SXY,10), colorbar, title('stress sxy')

% Principal stress magnitudes, maximum shear stress
S1 = 0.5*(SXX+SYY) + sqrt(0.25*(SXX-SYY).^2 + SXY.^2);
S2 = 0.5*(SXX+SYY) - sqrt(0.25*(SXX-SYY).^2 + SXY.^2);
SS = 0.5*(S1 - S2);
SM = 0.5*(S1 + S2);
figure, subplot(2,2,1), contourf(X,Y,S1,10), colorbar, title('stress s1')
subplot(2,2,2), contourf(X,Y,S2,10), colorbar, title('stress s2')
subplot(2,2,3), contourf(X,Y,SS,10), colorbar, title('stress ss')
subplot(2,2,4), contourf(X,Y,SM,10), colorbar, title('stress sm')

% Principal stress trajectories
x = linspace(0,2.5,26)+eps; y = linspace(0,2.5,26);
[X,Y] = meshgrid(x,y);
[TH,R] = cart2pol(X,Y);

ST = sin(TH); S2T = sin(2*TH); ST2 = ST.^2; 
CT = cos(TH); C2T = cos(2*TH); CT2 = CT.^2;
R2 = (ri./R).^2; R4 = R2.^2;

% Polar stress components
SRR = -(0.5*(sH+sh)*(1-R2))-(pm*R2)-(0.5*(sH-sh)*((1-4*R2+3*R4).*C2T));
STT = -(0.5*(sH+sh)*(1+R2))+(pm*R2)+(0.5*(sH-sh)*((1+3*R4).*C2T));
SRT = 0.5*(sH-sh)*((1+2*R2-3*R4).*S2T);
SRR(find(R<1)) = nan; STT(find(R<1)) = nan; SRT(find(R<1)) = nan;

% Cartesian stress components
SXX = SRR.*CT2+STT.*ST2-2*SRT.*CT.*ST;
SYY = SRR.*ST2+STT.*CT2+2*SRT.*CT.*ST;
SXY = (SRR-STT).*CT.*ST+SRT.*(CT2-ST2);

% Principal stress trajectories
G1 = 0.5*atan2(2*SXY, SXX-SYY); G1(find(R<1))=nan;
U1 = cos(G1); V1 = sin(G1);
quiver(X,Y,U1,V1,0.2,'.'), axis equal, hold on
G1 = G1+pi/2; U1 = cos(G1); V1 = sin(G1);
quiver(X,Y,U1,V1,0.4,'.')
G1 = G1+pi/2; U1 = cos(G1); V1 = sin(G1);
quiver(X,Y,U1,V1,0.2,'.')
G1 = G1+pi/2; U1 = cos(G1); V1 = sin(G1);
quiver(X,Y,U1,V1,0.4,'.')