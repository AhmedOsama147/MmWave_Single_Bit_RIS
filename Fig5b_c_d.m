% Parameters
close all; clc; clear;
M = [10 20 30]; N = [10 20 30]; % 32x32 RIS
theta_i = 0; phi_i = 0; % Broadside illumination
theta_d = 1:0.5:60; phi_d = 0; % Desired reflection angle
error_d1 = zeros(size(theta_d));
error_d2 = zeros(size(theta_d));
theta_HPBW2 = zeros(size(theta_d));
theta_HPBW1 = zeros(size(theta_d));
lambda = 1; % Normalized wavelength
k0 = 2*pi/lambda;

% Grid setup
dx = 0.5; dy = 0.5; % Element spacing (lambda/2)

% Phase calculation (traditional RIS)
for m = 1:length(M)
[X, Y] = meshgrid(-M(m)/2:M(m)/2-1, -N(m)/2:N(m)/2-1);

    for i = 1:length(theta_d)
    phi_0 = -k0 * (X*dx*sind(theta_d(i))*cosd(phi_d) + Y*dy*sind(theta_d(i))*sind(phi_d));
    phi_0 = wrapTo180(phi_0*(180/(pi)));
    
    phi_rand = (rand(M(m),N(m)))*180;
    
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
    
    for x = 1:M(m)
        for y = 1:N(m)
           AF1 = AF1 + exp(-1j*phi_total(x,y)*pi/180) * exp(-1j*k0*((x-M(m)/2-1)*dx*u + (y-N(m)/2-1)*dy*v));
           AF2 = AF2 + exp(-1j*phi_quant2(x,y)*pi/180) * exp(-1j*k0*((x-M(m)/2-1)*dx*u + (y-N(m)/2-1)*dy*v));
        end
    end
    
    AF1 = abs(AF1) / max(abs(AF1(:))); % Normalize
    AF2 = abs(AF2) / max(abs(AF2(:))); % Normalize
    
    [~,IND1] = max(AF1); ex_theta_d1 = abs(theta(IND1));
    [~,IND2] = max(AF2); ex_theta_d2 = abs(theta(IND2));   
    error_d1(i) = abs(theta_d(i)-ex_theta_d1);
    error_d2(i) = abs(theta_d(i)-ex_theta_d2);
    theta_HPBW1(i) = HPBW(theta,10*log10(AF1));
    theta_HPBW2(i) = HPBW(theta,10*log10(AF2));

end
subplot(2,2,m)
stem(theta_d,error_d1./theta_HPBW1,'LineWidth',0.75,'Color','r','Marker','none');
hold on
stem(theta_d,error_d2./theta_HPBW2,'LineWidth',0.5,'Color','k','Marker','none'); 
legend('$P^3$ RIS','Traditional RIS','Interpreter','Latex');
grid on;
xlabel('Reflective Angle');
ylabel('Normalized BPE');
Title = int2str(m*10) + "x" + int2str(m*10);
title(Title);
hold off
end