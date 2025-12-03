% Test Points 8 and 9 with MATLAB TDstressHS (half-space)
clear all;

% Triangle vertices
P1 = [-1.0, -1.0, -5.0];
P2 = [1.0, -1.0, -5.0];
P3 = [-1.0, 1.0, -4.0];

% Slip components
Ss = 1.0;   % Strike-slip
Ds = -1.0;  % Dip-slip
Ts = 2.0;   % Tensile-slip

% Elastic parameters
mu = 3.0e10;
lambda = 3.0e10;

fprintf('Testing Points 8 and 9 with MATLAB TDstressHS (Half-Space)\n');
fprintf('============================================================\n\n');

% Point 8
X = 3.0;
Y = -3.0;
Z = -6.0;
[Stress, Strain] = TDstressHS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Point 8: (%.1f, %.1f, %.1f)\n', X, Y, Z);
fprintf('  Exx = %.15e\n', Strain(1));
fprintf('  Eyy = %.15e\n', Strain(2));
fprintf('  Ezz = %.15e\n', Strain(3));
fprintf('\n');

% Point 9
X = -3.0;
Y = 3.0;
Z = -3.0;
[Stress, Strain] = TDstressHS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Point 9: (%.1f, %.1f, %.1f)\n', X, Y, Z);
fprintf('  Exx = %.15e\n', Strain(1));
fprintf('  Eyy = %.15e\n', Strain(2));
fprintf('  Ezz = %.15e\n', Strain(3));
fprintf('\n');
