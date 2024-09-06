% IEEE Transactions on Industrial Cyber Physical System
% WheelModel.m
% Created by Muhammad Hilmi
% Last Update: 06-09-2024

function [xw, yw] = WheelModel(x, y, yaw)

    % Calculate the four corner points of the wheel based on the yaw angle
    x1 = x - 0.6*cos(yaw) + 0.4*sin(yaw);
    x2 = x + 0.6*cos(yaw) + 0.4*sin(yaw);
    x3 = x + 0.6*cos(yaw) - 0.4*sin(yaw);
    x4 = x - 0.6*cos(yaw) - 0.4*sin(yaw);

    y1 = y - 0.6*sin(yaw) - 0.4*cos(yaw);
    y2 = y + 0.6*sin(yaw) - 0.4*cos(yaw);
    y3 = y + 0.6*sin(yaw) + 0.4*cos(yaw);
    y4 = y - 0.6*sin(yaw) + 0.4*cos(yaw);

    % Combine the points into vectors for output
    xw = [x1, x2, x3, x4];
    yw = [y1, y2, y3, y4];

end
