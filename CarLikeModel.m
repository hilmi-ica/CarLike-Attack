% IEEE Transactions on Industrial Cyber Physical System
% CarLikeModel.m
% Created by Muhammad Hilmi
% Last Update: 06-09-2024

function H = CarLikeModel(x, y, teta, delta, H)

    % Adjust delta for wheel angle
    delta = teta + delta;

    % Define the vehicle's corner positions
    xCorners = [x - 1.5*cos(teta) + 1*sin(teta), x - 1*cos(teta) + 1.5*sin(teta), ...
                x + 4*cos(teta) + 1.5*sin(teta), x + 4.5*cos(teta) + 1*sin(teta), ...
                x + 4.5*cos(teta) - 1*sin(teta), x + 4*cos(teta) - 1.5*sin(teta), ...
                x - 1*cos(teta) - 1.5*sin(teta), x - 1.5*cos(teta) - 1*sin(teta)];
    
    yCorners = [y - 1.5*sin(teta) - 1*cos(teta), y - 1*sin(teta) - 1.5*cos(teta), ...
                y + 4*sin(teta) - 1.5*cos(teta), y + 4.5*sin(teta) - 1*cos(teta), ...
                y + 4.5*sin(teta) + 1*cos(teta), y + 4*sin(teta) + 1.5*cos(teta), ...
                y - 1*sin(teta) + 1.5*cos(teta), y - 1.5*sin(teta) + 1*cos(teta)];

    % Calculate front and rear wheel positions
    xpflw = x + 3*cos(teta) - 1.5*sin(teta);
    ypflw = y + 3*sin(teta) + 1.5*cos(teta);
    xpfrw = x + 3*cos(teta) + 1.5*sin(teta);
    ypfrw = y + 3*sin(teta) - 1.5*cos(teta);

    [xflw, yflw] = WheelModel(xpflw, ypflw, delta);
    [xfrw, yfrw] = WheelModel(xpfrw, ypfrw, delta);

    xprlw = x + 0*cos(teta) - 1.5*sin(teta);
    yprlw = y + 0*sin(teta) + 1.5*cos(teta);
    xprrw = x + 0*cos(teta) + 1.5*sin(teta);
    yprrw = y + 0*sin(teta) - 1.5*cos(teta);

    [xrlw, yrlw] = WheelModel(xprlw, yprlw, teta);
    [xrrw, yrrw] = WheelModel(xprrw, yprrw, teta);

    % Plot the car's body
    for i = 1:8
        j = mod(i, 8) + 1; % wrap around to connect the last point to the first
        H(i) = plot([xCorners(i) xCorners(j)], [yCorners(i) yCorners(j)], 'k-', 'LineWidth', 2);
    end

    % Plot the wheels
    wheelPairs = {xflw, yflw; xfrw, yfrw; xrlw, yrlw; xrrw, yrrw};
    offset = 8;
    for k = 1:4
        [xw, yw] = wheelPairs{k, :};
        for i = 1:4
            j = mod(i, 4) + 1; % wrap around to connect the last point to the first
            H(offset + (k-1)*4 + i) = plot([xw(i) xw(j)], [yw(i) yw(j)], '-', 'Color', "#0072BD", 'LineWidth', 2);
        end
    end
end
