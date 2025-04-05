function [labi] = lambda_i(V, rotor_speed, R, rho, S, Cdf, W)
    mu = V / (rotor_speed * R);
    D = 0.5 * rho * V^2 * S * Cdf; % Drag force (N)
    
    CT1 = sqrt(W^2 + D^2) / (rho * (rotor_speed * R)^2 * pi * R^2);
    
    % Solve for lambda_i using fzero
    lambda_i_guess = 0.05;
    func = @(lambda_i) CT1 - 2 * lambda_i * sqrt((V / (rotor_speed * R))^2 + (V * (D/W) / (rotor_speed * R) + lambda_i)^2);
    labi = fzero(func, lambda_i_guess);
end