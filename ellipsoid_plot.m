function [x_rot, y_rot, z_rot] = ellipsoid_plot(P,x0)
% Eigenvalues and eigenvectors
[V, D] = eig(P);

% Semi-axes lengths based on the eigenvalues
a = 1 / sqrt(D(1,1));
b = 1 / sqrt(D(2,2));
c = 1 / sqrt(D(3,3));

% Generate an axis-aligned ellipsoid using the built-in function
[x, y, z] = ellipsoid(0, 0, 0, a, b, c, 30);

% Apply the rotation using the eigenvectors
ellipsoid_coords = V * [x(:)'; y(:)'; z(:)'];
x_rot = reshape(ellipsoid_coords(1, :), size(x));
y_rot = reshape(ellipsoid_coords(2, :), size(y));
z_rot = reshape(ellipsoid_coords(3, :), size(z));

% Plot the rotated ellipsoid
surf(x_rot+x0(1), y_rot+x0(2), z_rot+x0(3), 'FaceAlpha', 0.6,'EdgeAlpha', 0.3);
hold on;
grid on;
end