% IEEE Transactions on Industrial Cyber Physical System
% CarLikeSim.m
% Created by Muhammad Hilmi
% Last Update: 18-10-2025

clear all;
close all;

% Load Data
load("Odometry.mat");
load("Noise.mat");
load("RBFNN.mat");

% Initialize Parameters
l = 3; dt = 0.05; vm = 2;
xi = 50; yi = 0; tetai = pi/2;

% System Matrices
A = zeros(5); B = [0 0; 0 0; 0 0; 1 0; 0 1]; C = [1 0 0 0 0; 0 1 0 0 0];
X = eye(2); Af = -10*eye(2);
Aa = [A zeros(5, 2); -Af*C Af]; Ba = [B; zeros(2)]; Ca = [zeros(2, 5) eye(2)];
Xa = [zeros(5, 2); -Af*X]; Aa = eye(7) + Aa*dt;

% Initial States
x = [xi yi tetai vm 0 0 0]';
xnoise = x; xnlo = x; xhat = x; xrbf = x; xreal = [xi yi tetai vm 0]';
xArray = []; xnoiseArray = []; xnloArray = []; xhatArray = []; xrbfArray = []; xrealArray = [];
yreal = [xi yi]'; yrealArray = [];

theta = [0; 0]; thetanlo = [0; 0]; thetahat = [0; 0]; thetarbf = [0; 0];
thetaArray = []; thetanloArray = []; thetahatArray = []; thetarbfArray = [];

% Parameters for Estimation
v = 0.01; e1 = 1; e2 = 100; p = -10; b = 0.01; n = e1/e2;
a1 = v + e1*p + e2*b; a2 = (n*e2 - e1)/2; a3 = a2; a4 = -e2;
J = 3*eye(7); Y = [0 0 0 0 0 1 0; 0 0 0 0 0 0 1]; Q = 0.01*eye(7);

eivals = eig([-J + a1*eye(7), Aa'*J - Ca'*Y + a2*eye(7); J*Aa - Y'*Ca + a3*eye(7), J + a4*eye(7)]);
alpha = 0.95; lambda = 0.999; Pk = 0.01*eye(7); Qk = 0.01*eye(7); Rk = 10^8*eye(2);
Kappa = zeros(7, 2); Sk = 0.01*eye(2);

Pk2 = 0.01*eye(7); Qk2 = 0.01*eye(7); Rk2 = 10^8*eye(2);
Kappa2 = zeros(7, 2); Sk2 = 0.01*eye(2);

x_prev = x; xnoise_prev = x + [zeros(5,1); -Af*zeros(2,1)];
xnlo_prev = xnlo; thetahat_prev = thetahat;
at = 62;

% Main Simulation Loop
for i = 1:length(a)
    % Attack Scenarios
    if i*dt > 0
        theta(1) = 50 * max(-1, min(0, (at - i*dt)));
        theta(2) = 20 * exp(-(i*dt - 65)^2 / 2*0.05^2);
    end
    if i*dt > 1240
        theta(1) = ((i - 2500) * 0.0005) * 10;
    end
    u = [a(i); w(i)];

    % Store Data for Plotting
    xArray = [xArray, x]; xnoiseArray = [xnoiseArray, xnoise];
    xnloArray = [xnloArray, xnlo]; xhatArray = [xhatArray, xhat];
    xrbfArray = [xrbfArray, xrbf]; xrealArray = [xrealArray, xreal];
    yrealArray = [yrealArray, yreal];
    thetaArray = [thetaArray, theta]; thetanloArray = [thetanloArray, thetanlo];
    thetahatArray = [thetahatArray, thetahat]; thetarbfArray = [thetarbfArray, thetarbf];

    % System Dynamics and Estimation Updates
    fx = [x(4)*cos(x(3)); x(4)*sin(x(3)); x(4)/l*tan(x(5)); 0; 0; 0; 0];
    xnoise = Aa*x + fx*dt + Ba*u*dt + Xa*theta*dt + dArray(:,i);
    ynoise = Ca*xnoise;
    x = Aa*x + fx*dt + Ba*u*dt + Xa*theta*dt;
    y = Ca*x;

    % Real Dynamics
    fxreal = [xreal(4)*cos(xreal(3)); xreal(4)*sin(xreal(3)); xreal(4)/l*tan(xreal(5)); 0; 0];
    xreal = eye(5)*xreal + fxreal*dt + B*u*dt;
    yreal = C*xreal + X*theta + dArray(6:7,i)/10;

    % RBFNN Estimation
    trbf = tic;
    results = sim(rbf, [u; xnoise]);
    xrbf = results(1:7,:);
    thetarbf = results(8:9,:);
    delrbf(i) = toc(trbf);

    % AXKF Estimation
    txkf = tic;
    K = Aa - inv(J)*Y'*Ca;
    P = 2*inv(Xa*Xa' + Q);
    fxnlo = [xnlo(4)*cos(xnlo(3)); xnlo(4)*sin(xnlo(3)); xnlo(4)/l*tan(xnlo(5)); 0; 0; 0; 0];
    xnlo = Aa*xnlo + fxnlo*dt + Ba*u*dt + Xa*thetanlo*dt + inv(J)*Y'*(ynoise - Ca*xnlo)*dt;
    thetanlo = thetanlo + Xa'*P*((xnoise - xnlo) - K*(xnoise_prev - xnlo_prev) - (fx - fxnlo));
    Fxh = [zeros(5,5) zeros(5,2); zeros(2,7)];
    Fxh(1:5, 1:5) = [0 0 -xnlo(4)*sin(xnlo(3)) cos(xnlo(3)) 0;
                     0 0 xnlo(4)*cos(xnlo(3)) sin(xnlo(3)) 0;
                     0 0 0 tan(xnlo(5))/l (xnlo(4)/l)*(sec(xnlo(5)))^2;
                     0 0 0 0 0;
                     0 0 0 0 0];
    Pk = (Fxh*dt)*Pk*(Fxh*dt)' + Qk;
    Sigma = Ca*Pk*Ca' + Rk;
    Kk = Pk*Ca'*inv(Sigma);
    Pk = (eye(7) - Kk*Ca)*Pk;
    Omega = Ca*(Fxh*dt)*Kappa + Ca*Xa*dt;
    Kappa = (eye(7) - Kk*Ca)*(Fxh*dt)*Kappa + (eye(7) - Kk*Ca)*Xa*dt;
    Lambda = inv(lambda*Sigma + Omega*Sk*Omega');
    Gamma = Sk*Omega'*Lambda;
    Sk = Sk/lambda - Sk*Omega'*Lambda*Omega*Sk/lambda;
    ytilde = ynoise - Ca*xhat;
    thetahat = thetanlo + Gamma*ytilde;
    xhat = Aa*xnlo + Fxh*xhat*dt + fxnlo*dt - Fxh*xnlo*dt + Ba*u*dt + Xa*thetahat_prev*dt + Kk*ytilde*dt + Kappa*(thetahat - thetahat_prev)*dt;
    delxkf(i) = toc(txkf);

    % Update Previous States
    x_prev = x; xnoise_prev = xnoise; xnlo_prev = xnlo; thetahat_prev = thetahat;

    % Error Calculations
    exhat(:,i) = xnoise - xhat; 
    exrbf(:,i) = xnoise - xrbf;
end

% Define time vector for plotting
t = 0:dt:i*dt;
t = t(1:end-1);

% Plot theta estimation comparison
fh = figure(1);
fh.Position = [0 50 1000 450];

% Plot for theta_1
subplot(2, 1, 1)
plot(t, thetaArray(1,:), 'k', 'LineWidth', 3)
hold on
plot(t, thetarbfArray(1,:), '--', 'Color', "#EDB120", 'LineWidth', 3)
plot(t, thetahatArray(1,:), '--', 'Color', "#77AC30", 'LineWidth', 3)
legend('\theta_{1,Actual}', '\theta_{1,RBFNN}', '\theta_{1,AXKF}', 'Location', 'northwest')
grid on; grid minor
xlim('tight'); ylim('tight')
ylabel('\theta_1 (m)')
set(gca, 'Color', 'white', 'FontSize', 14)
hold off

% Plot for theta_2
subplot(2, 1, 2)
plot(t, thetaArray(2,:), 'k', 'LineWidth', 3)
hold on
plot(t, thetarbfArray(2,:), '--', 'Color', "#EDB120", 'LineWidth', 3)
plot(t, thetahatArray(2,:), '--', 'Color', "#77AC30", 'LineWidth', 3)
legend('\theta_{2,Actual}', '\theta_{2,RBFNN}', '\theta_{2,AXKF}', 'Location', 'northwest')
grid on; grid minor
xlabel('Time (s)')
xlim('tight'); ylim('tight')
ylabel('\theta_2 (m)')
set(gcf, 'Color', 'white')
set(gca, 'FontSize', 14)
hold off

% Initialize Car-Like Model plot
H = zeros(1,24);
n = length(xArray(1,:));
fh = figure(2);
fh.Position = [0 50 1000 800];

% Plot state estimation comparison
plot(xnoiseArray(1, :), xnoiseArray(2, :), 'k-', 'LineWidth', 3)
hold on
plot(yrealArray(1, :), yrealArray(2, :), ':', 'Color', "#A2142F", 'LineWidth', 3)
plot(xrbfArray(1, :), xrbfArray(2, :), '--', 'Color', "#EDB120", 'LineWidth', 3)
plot(xhatArray(1, :), xhatArray(2, :), '--', 'Color', "#77AC30", 'LineWidth', 3)

% Add car-like model markers at start, middle, and end
H = CarLikeModel(xArray(1, 1), xArray(2, 1), xArray(3, 1), xArray(5, 1));
H = CarLikeModel(xArray(1, ceil(n/2)), xArray(2, ceil(n/2)), xArray(3, ceil(n/2)), xArray(5, ceil(n/2)));
H = CarLikeModel(xArray(1, n), xArray(2, n), xArray(3, n), xArray(5, n));

% Add start and stop text annotations
text(xnoiseArray(1, 1) - 2, xnoiseArray(2, 1) - 3, 'START', 'Color', 'k', 'FontSize', 10); 
text(xnoiseArray(1, end) - 0.5, xnoiseArray(2, end) - 3.5, 'STOP', 'Color', 'k', 'FontSize', 10);

% Add legend and axis formatting
legend('True State', 'Spoofed Measurement', 'RBFNN Estimation', 'AXKF Estimation', 'Location', 'northeast')
grid on; grid minor
axis([-57 57 -10 80])
xlabel('x (m)')
ylabel('y (m)')
set(gcf, 'Color', 'white')
set(gca, 'FontSize', 14)
hold off
