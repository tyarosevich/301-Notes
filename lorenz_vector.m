function dy = lorenz_vector(t, y, sigma, beta, rho)
% y is an indexed set of meshgrid values
dy = [
    sigma * (y(2,:,:,:) - y(1,:,:,:)); % sigma(y - x)
    y(1,:,:,:) .* (rho - y(3,:,:,:)) - y(2,:,:,:); % x*(rho - z)
    y(1,:,:,:) .* y(2,:,:,:) - beta * y(3,:,:,:); % x * y - beta * z
    ];
%this is just vectorized