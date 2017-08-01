%% Coding up second order accurate difference.

function [x_out, dy_dx] = backdiff_2nd(f, a, b, N)

x = linspace(a,b,N);
dx = (b-a)/(N-1);
y = f(x);
dy_dx = (3*y(3:end) - 4*y(2:end-1) + y(1:end-2))/(2*dx);
x_out = x(1:end-1);
end 