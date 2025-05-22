% Parameters
close all; clc; clear;
M = 10; N = 10; % 10x10 RIS
theta_i = 0; phi_i = 0; % Broadside illumination
theta_d = 4; phi_d = 0; % Desired reflection angle
lambda = 1; % Normalized wavelength
k0 = 2*pi/lambda;

% Grid setup
[x, y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
dx = 0.5; dy = 0.5; % Element spacing (lambda/2)

% Phase calculation (traditional RIS)
phi_0 = -k0 * (x*dx*sind(theta_d)*cosd(phi_d) + y*dy*sind(theta_d)*sind(phi_d));
phi_0 = wrapTo180(phi_0*(180/(pi)));

phi_rand = (rand(M,N))*180;

phi_mn = phi_0 - phi_rand ;
phi_quant1 = abs(180 * floor((phi_mn)/180 + 0.5)) ; % 1-bit quantization
phi_quant2 = abs(180 * floor((phi_0)/180 + 0.5)) ; % 1-bit quantization
phi_total = wrapTo180(phi_quant1 + phi_rand);
% Radiation pattern (array factor)
% Angular grid
theta = linspace(-90, 90, 512);
phi = linspace(-180, 180, 512);
[Theta, Phi] = meshgrid(theta, phi);
u = sind(theta).*cosd(90);
v = sind(theta).*sind(90);
AF1 = 0;
AF2 = 0;

for x = 1:M
    for y = 1:N
       AF1 = AF1 + exp(-1j*phi_total(x,y)*pi/180) * exp(-1j*k0*((x-M/2-1)*dx*u + (y-N/2-1)*dy*v));
       AF2 = AF2 + exp(-1j*phi_quant2(x,y)*pi/180) * exp(-1j*k0*((x-M/2-1)*dx*u + (y-N/2-1)*dy*v));
    end
end

AF1 = abs(AF1) / max(abs(AF1(:))); % Normalize
AF2 = abs(AF2) / max(abs(AF2(:))); % Normalize

% Plot
figure;
plot(theta, 10*log10(AF1),'Color','r','LineWidth',2);
hold on
plot(theta, 10*log10(AF2),'Color','k','LineWidth',2); 
view(0,90);
xlabel('Reflection Angle in degrees');
ylabel('Normalized Magnitude in (dB)')
xlim([-90 90]);
ylim([-30 0]);
xticks(-90:30:90);
yticks(-30:10:0);
legend('P^3 RIS','Traditional RIS')