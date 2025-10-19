% IEEE Transactions on Industrial Cyber Physical System
% CarLikeMC.m
% Created by Muhammad Hilmi
% Last Update: 20-10-2025

clear; close all; clc;

% Load Data
load("Odometry.mat");
load("Noise.mat");

% Fixed Parameters
l = 3; dt = 0.05; vm = 2;
xi = 50; yi = 0; tetai = pi/2;

A = zeros(5); 
B = [0 0; 0 0; 0 0; 1 0; 0 1]; 
C = [1 0 0 0 0; 0 1 0 0 0];
X = eye(2); Af = -10*eye(2);
Aa = [A zeros(5,2); -Af*C Af]; 
Ba = [B; zeros(2)]; 
Ca = [zeros(2,5) eye(2)];
Xa = [zeros(5,2); -Af*X]; 
Aa = eye(7) + Aa*dt;

% Constants for LMI
v = 0.01; e1 = 1; e2 = 100; p = -10; b = 0.01; n = e1/e2;
a1 = v + e1*p + e2*b; 
a2 = (n*e2 - e1)/2; 
a3 = a2; 
a4 = -e2;

H0 = 3*eye(7); J0 = 0.01*eye(7);
Y = [0 0 0 0 0 1 0; 0 0 0 0 0 0 1];

% Monte Carlo Grid
kH_vals = linspace(0.1, 10, 30);
kJ_vals = linspace(0.1, 10, 30);

nxhat_map = nan(length(kH_vals),length(kJ_vals));
nthat_map = nan(length(kH_vals),length(kJ_vals));
feasible = false(size(nxhat_map));

% Monte Carlo Sweep
for ii = 1:length(kH_vals)
    for jj = 1:length(kJ_vals)
        kH = kH_vals(ii);
        kJ = kJ_vals(jj);

        % Construct SPD-like H, J
        H = abs(kH)*H0;
        J = abs(kJ)*J0;

        % LMI feasibility check
        LHS = [ -H + a1*eye(7), H + a2*eye(7), H;
                 H + a3*eye(7),  H + a4*eye(7), zeros(7);
                 H, zeros(7), -H ];
        LHS = (LHS + LHS')/2;
        eigL = eig(LHS);

        if all(eigL < -1e-6)
            feasible(ii,jj) = true;
            
            % Initialize States
            x = [xi yi tetai vm 0 0 0]';
            xnoise = x; xnlo = x; xhat = x;
            xreal = [xi yi tetai vm 0]';
            theta = [0; 0]; thetanlo = [0; 0]; thetahat = [0; 0];

            Pk = 0.01*eye(7); Qk = 0.01*eye(7); Rk = 10^8*eye(2);
            Upsilon = zeros(7,2); Sk = 0.01*eye(2);
            x_prev = x; xnoise_prev = x; xnlo_prev = xnlo; thetahat_prev = thetahat;
            at = 62;

            for i = 1:length(a)
                % Attack
                if i*dt > 0
                    theta(1) = 50 * max(-1, min(0, (at - i*dt)));
                    theta(2) = 20 * exp(-(i*dt - 65)^2 / 2*0.05^2);
                end
                if i*dt > 1240
                    theta(1) = ((i - 2500) * 0.0005) * 10;
                end
                u = [a(i); w(i)];

                % Dynamics + Noise
                fx = [x(4)*cos(x(3)); x(4)*sin(x(3)); x(4)/l*tan(x(5)); 0; 0; 0; 0];
                xnoise = Aa*x + fx*dt + Ba*u*dt + Xa*theta*dt + dArray(:,i);
                ynoise = Ca*xnoise;
                x = Aa*x + fx*dt + Ba*u*dt + Xa*theta*dt;

                % Real dynamics
                fxr = [xreal(4)*cos(xreal(3)); xreal(4)*sin(xreal(3)); xreal(4)/l*tan(xreal(5)); 0; 0];
                xreal = eye(5)*xreal + fxr*dt + B*u*dt;
                y = Ca*x;

                % AXKF Estimation
                K = Aa - inv(H)*Y'*Ca;
                P = 2*inv(Xa*Xa' + J);
                fxnlo = [xnlo(4)*cos(xnlo(3)); xnlo(4)*sin(xnlo(3)); xnlo(4)/l*tan(xnlo(5)); 0; 0; 0; 0];
                xnlo = Aa*xnlo + fxnlo*dt + Ba*u*dt + Xa*thetanlo*dt + inv(H)*Y'*(ynoise - Ca*xnlo)*dt;
                thetanlo = thetanlo + Xa'*P*((xnoise - xnlo) - K*(xnoise_prev - xnlo_prev) - (fx - fxnlo));

                Fxh = [zeros(5,5) zeros(5,2); zeros(2,7)];
                Fxh(1:5, 1:5) = [0 0 -xnlo(4)*sin(xnlo(3)) cos(xnlo(3)) 0;
                                 0 0 xnlo(4)*cos(xnlo(3)) sin(xnlo(3)) 0;
                                 0 0 0 tan(xnlo(5))/l (xnlo(4)/l)*(sec(xnlo(5)))^2;
                                 0 0 0 0 0;
                                 0 0 0 0 0];
                Pk = (Fxh*dt)*Pk*(Fxh*dt)' + Qk;
                Sigma = Ca*Pk*Ca' + Rk;
                Kk = Pk*Ca'/Sigma;
                Pk = (eye(7) - Kk*Ca)*Pk;
                Omega = Ca*(Fxh*dt)*Upsilon + Ca*Xa*dt;
                Upsilon = (eye(7) - Kk*Ca)*(Fxh*dt)*Upsilon + (eye(7) - Kk*Ca)*Xa*dt;
                Lambda = inv(0.999*Sigma + Omega*Sk*Omega');
                Gamma = Sk*Omega'*Lambda;
                Sk = Sk/0.999 - Sk*Omega'*Lambda*Omega*Sk/0.999;
                ytilde = ynoise - Ca*xhat;
                thetahat = thetanlo + Gamma*ytilde;
                xhat = Aa*xnlo + Fxh*xhat*dt + fxnlo*dt - Fxh*xnlo*dt + ...
                       Ba*u*dt + Xa*thetahat_prev*dt + Kk*ytilde*dt + Upsilon*(thetahat - thetahat_prev)*dt;

                % Update
                x_prev = x; xnoise_prev = xnoise; xnlo_prev = xnlo; thetahat_prev = thetahat;
                exhat(:,i) = xnoise - xhat; 
                ethat(:,i) = theta - thetahat;
            end

            nxhat_map(ii,jj) = norm(sqrt(mean(exhat.^2,2)));
            nthat_map(ii,jj) = norm(sqrt(mean(ethat.^2,2)));

        else
            nxhat_map(ii,jj) = NaN; 
            nthat_map(ii,jj) = NaN;
        end
    end
end

% Visualization
close all
[KH,KJ] = meshgrid(kH_vals,kJ_vals);
Z = nxhat_map'; C = nthat_map';

figure('Color','w','Position',[100 100 900 700]);
surf(KH, KJ, Z, C, 'LineWidth', 2);
hold on; cb = colorbar('southoutside','LineWidth',1.5); view(-30,30);
xlabel('k_1'); ylabel('k_2'); zlabel('||{\bfz\rm}_{RMSE}||');
ylabel(cb, '||{\bf\theta\rm}_{RMSE}||', 'FontSize', 14)

% Mark unfeasible regions
unfidx = find(~feasible(:,1)');
kHuf = [kH_vals(unfidx(1)) kH_vals(unfidx(end))];
kJuf = [kJ_vals(1) kJ_vals(end)];
nXuf = [min(nxhat_map,[],"all") max(nxhat_map,[],"all")];
Xunf = [kHuf(1) kHuf(1) kHuf(1) kHuf(2) kHuf(2) kHuf(2); ...
    kHuf(1) kHuf(2) kHuf(1) kHuf(2) kHuf(1) kHuf(1); ...
    kHuf(1) kHuf(2) kHuf(2) kHuf(2) kHuf(1) kHuf(1); ...
    kHuf(1) kHuf(1) kHuf(2) kHuf(2) kHuf(2) kHuf(2)];
Yunf = [kJuf(1) kJuf(1) kJuf(1) kJuf(1) kJuf(2) kJuf(2); ...
    kJuf(2) kJuf(1) kJuf(2) kJuf(2) kJuf(2) kJuf(2); ...
    kJuf(2) kJuf(1) kJuf(2) kJuf(2) kJuf(2) kJuf(1); ...
    kJuf(1) kJuf(1) kJuf(1) kJuf(1) kJuf(2) kJuf(1)];
Zunf = [nXuf(1) nXuf(1) nXuf(1) nXuf(1) nXuf(1) nXuf(2); ...
    nXuf(1) nXuf(1) nXuf(1) nXuf(1) nXuf(1) nXuf(2); ...
    nXuf(2) nXuf(2) nXuf(1) nXuf(2) nXuf(2) nXuf(2); ...
    nXuf(2) nXuf(2) nXuf(1) nXuf(2) nXuf(2) nXuf(2)];
fill3(Xunf,Yunf,Zunf,'red', LineWidth=1.5)

legend({'','Infeasible region'}); grid on; box on; hold off;
ax = gca; ax.LineWidth = 1.5; ax.GridLineWidth = 1;
ax.FontSize = 14; set(gcf, 'Color', 'white')