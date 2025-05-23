% Parameters
close all; clc; clear;
M = 10; N = 10; % 10x10 RIS
theta_i = 0; phi_i = 0; % Broadside illumination
theta_d = 1:60; phi_d = 0; % Desired reflection angle
ex_theta_d1 = zeros(size(theta_d));
ex_theta_d2 = zeros(size(theta_d));
lambda = 1; % Normalized wavelength
k0 = 2*pi/lambda;

% Grid setup
[X, Y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
dx = 0.5; dy = 0.5; % Element spacing (lambda/2)

% Phase calculation (traditional RIS)
for i = theta_d
    phi_0 = -k0 * (X*dx*sind(i)*cosd(phi_d) + Y*dy*sind(i)*sind(phi_d));
    phi_0 = wrapTo180(phi_0*(180/(pi)));
    
    phi_rand = (rand(M,N))*180;
    
    phi_mn = phi_0 - phi_rand ;
    phi_quant1 = abs(180 * floor((phi_mn)/180 + 0.5)) ; % 1-bit quantization
    phi_quant2 = abs(180 * floor((phi_0)/180 + 0.5)) ; % 1-bit quantization
    phi_total = wrapTo180(phi_quant1 + phi_rand);
    % Radiation pattern (array factor)
    % Angular grid
    theta = linspace(-90, 90, 512);
    u = sind(theta).*cosd(90);
    v = sind(theta).*sind(90);
    AF1 = 0; % AF for P3-RIS
    AF2 = 0; % AF for traditional-RIS
    
    for x = 1:M
        for y = 1:N
           AF1 = AF1 + exp(-1j*phi_total(x,y)*pi/180) * exp(-1j*k0*((x-M/2-1)*dx*u + (y-N/2-1)*dy*v));
           AF2 = AF2 + exp(-1j*phi_quant2(x,y)*pi/180) * exp(-1j*k0*((x-M/2-1)*dx*u + (y-N/2-1)*dy*v));
        end
    end
    
    AF1 = abs(AF1) / max(abs(AF1(:))); % Normalize
    AF2 = abs(AF2) / max(abs(AF2(:))); % Normalize
    
    [~,IND1] = max(AF1); ex_theta_d1(i) = abs(theta(IND1));
    [~,IND2] = max(AF2); ex_theta_d2(i) = abs(theta(IND2));   
end

plot(theta_d,ex_theta_d1,'LineWidth',2,'Color','r');
hold on
plot(theta_d,ex_theta_d2,'LineWidth',2,'Color','k');
hold off
legend('$P^3$ RIS','Traditional RIS','Interpreter','Latex');
grid on;
xlim([0 60]);
ylim([0 60]);
xlabel('Desired Beam Direction (deg)');
ylabel('Actual Beam Direction (deg)');