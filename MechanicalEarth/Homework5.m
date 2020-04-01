%Plot satilite data
D = load('Izmit_data.mat').Izmit_data;
% deepest extent of fault (km)
d_1 = 1;
% shallowest extent of fault (meters) (keep as 0)
d_2 = 0;
% magnitude of slip (meters)
s = 1;
%values for the x_1 axis (km)
x_1 = linspace(-80, 50, 50);
figure(1);

hold on
u_first = u_3(x_1, s, d_1, d_2);
plot(x_1, u_first);
d_1 = 10;
u_second = u_3(x_1, s, d_1, d_2);
plot(x_1, u_second);
d_1 = 20;
u_third = u_3(x_1, s, d_1, d_2);
plot(x_1, u_third);
plot(D(:,1), D(:,2));
title("Effects on Displacements of Varying Model Fault Depths(km)");
legend("1km", "10km", "20km", "Satillite");
xlabel("North-South km from fault");
ylabel("Displacement (m) in East-West direction");
hold off

figure(2);
hold on
d_1 = 10;
u_first = u_3(x_1, s, d_1, d_2);
plot(x_1, u_first);
s = 4;
u_second = u_3(x_1, s, d_1, d_2);
plot(x_1, u_second);
s = 7;
u_third = u_3(x_1, s, d_1, d_2);
plot(x_1, u_third);
plot(D(:,1), D(:,2));
legend("1m", "4m", "7m", "Satillite");
title("Effects on Displacements of Varying Model Fault Slip(m)");
xlabel("North-South km from fault");
ylabel("Displacement (m) in East-West direction");
hold off

figure(3);

s = 4;
d_1 =7;
u = u_3(x_1, s, d_1, d_2);
hold on
plot(x_1, u);
plot(D(:,1), D(:,2));
title("Comparison of Model with Slip(m) = 4 and Depth(km) = 7 to Satillite Data");
legend("Model", "Satillite");
xlabel("North-South km from fault");
ylabel("Displacement (m) in East-West direction");
hold off


%displacement along x2 = 0
function u = u_3(x_1, s, d_1, d_2)
u = -s/pi  * (atan(x_1/d_1) - atan(x_1/d_2));
end
