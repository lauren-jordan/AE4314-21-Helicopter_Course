function [a1, theta0] = compute_pitch(V, rotor_speed, R, rho, S, Cdf, W, cla, sigma)
    mu = V / (rotor_speed * R);
    D = 0.5 * rho * V^2 * S * Cdf; % Drag force (N)
    
    CT1 = sqrt(W^2 + D^2) / (rho * (rotor_speed * R)^2 * pi * R^2);
    
    % Solve for lambda_i using fzero
    lambda_i_guess = 0.05;
    func = @(lambda_i) CT1 - 2 * lambda_i * sqrt((V / (rotor_speed * R))^2 + (V * (D/W) / (rotor_speed * R) + lambda_i)^2);
    lambda_i = fzero(func, lambda_i_guess);
    
    % Solve system of equations
    A = [1 + (3/2) * mu^2, -(8/3) * mu;
         -mu, 1 + (2/3) * mu^2];
    
    B = [-2 * mu * lambda_i - 2 * mu^2 * (D/W);
          mu * (D/W) + lambda_i + 4 * CT1 / (cla * sigma)];
    
    X = A \ B; % Solve for a1 and theta0
    
    a1 = X(1); % Convert to degrees
    theta0 = X(2); % Convert to degrees
end