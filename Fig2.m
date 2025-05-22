% Parameters
close all; clc; clear;
M = 32; N = 32; % 32x32 RIS
theta_i = 0; phi_i = 0; % Broadside illumination
theta_d = 30; phi_d = 0; % Desired reflection angle
lambda = 1; % Normalized wavelength
k0 = 2*pi/lambda;

% Grid setup
[x, y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
dx = 0.5; dy = 0.5; % Element spacing (lambda/2)

% Phase calculation (traditional RIS)
phi_0 = -k0 * (x*dx*sind(theta_d)*cosd(phi_d) + y*dy*sind(theta_d)*sind(phi_d));
phi_0 = phi_0*(180/(pi));

phi_0 = wrapTo180(phi_0);
phi_quant = abs(180 * floor((phi_0)/180 + 0.5)) ; % 1-bit quantization

% Radiation pattern (array factor)
% Angular grid
theta = linspace(-90, 90, 512);
phi = linspace(-180, 180, 512);
[Theta, Phi] = meshgrid(theta, phi);
u = sind(Theta).*cosd(Phi);
v = sind(Theta).*sind(Phi);
AF = 0;
for x = 1:M
    for y = 1:N
        AF = AF + exp(-1j*phi_quant(x,y)*pi/180) * exp(-1j*k0*((x-M/2-1)*dx*u + (y-N/2-1)*dy*v));
    end
end
AF = abs(AF) / max(abs(AF(:))); % Normalize

% Plot
figure;
subplot(2,2,1); imagesc(x,y,zeros(M,N)); title('(a) \phi^{rand}_{MN} = 0'); 
axis square; colorbar; colormap jet; 
subplot(2,2,2); imagesc(phi_0); title('(b) excitation phase 32 × 32'); 
axis square; colorbar;
subplot(2,2,3); imagesc(phi_quant); title('(c) Quantized Phase \phi^{quant}_{mn}'); 
axis square; colorbar;
subplot(2,2,4); 
mesh(u, v, 10*log10(AF)); title('(d) normalized radiation pattern obtained in u − v plane'); 
axis equal;
c = colorbar;
colormap('jet');
view(0,90);