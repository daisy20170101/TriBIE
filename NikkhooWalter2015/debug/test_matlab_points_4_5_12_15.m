% Test Points 4, 5, 12, 15 with MATLAB TDstressHS
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

fprintf('============================================================\n');
fprintf('Testing Points 4, 5, 12, 15 with MATLAB TDstressHS\n');
fprintf('============================================================\n\n');

fprintf('Triangle vertices:\n');
fprintf('  P1 = (%.1f, %.1f, %.1f)\n', P1(1), P1(2), P1(3));
fprintf('  P2 = (%.1f, %.1f, %.1f)\n', P2(1), P2(2), P2(3));
fprintf('  P3 = (%.1f, %.1f, %.1f)\n', P3(1), P3(2), P3(3));
fprintf('\n');

% Point 4: (7.0, -1.0, -5.0)
fprintf('------------------------------------------------------------\n');
fprintf('Point 4: (7.0, -1.0, -5.0)\n');
fprintf('Same Y,Z as P2, but X=7.0 (6 units from P2 in X)\n');
[Stress, Strain] = TDstressHS(7.0, -1.0, -5.0, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('MATLAB Result: Exx = %.15e\n', Strain(1));
if isnan(Strain(1))
    fprintf('  -> NaN (singular)\n');
else
    fprintf('  -> Finite value\n');
end
fprintf('\n');

% Point 5: (-7.0, -1.0, -5.0)
fprintf('------------------------------------------------------------\n');
fprintf('Point 5: (-7.0, -1.0, -5.0)\n');
fprintf('Same Y,Z as P1, but X=-7.0 (6 units from P1 in X)\n');
[Stress, Strain] = TDstressHS(-7.0, -1.0, -5.0, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('MATLAB Result: Exx = %.15e\n', Strain(1));
if isnan(Strain(1))
    fprintf('  -> NaN (singular)\n');
else
    fprintf('  -> Finite value\n');
end
fprintf('\n');

% Point 12: (1.0, -1.0, -1.0)
fprintf('------------------------------------------------------------\n');
fprintf('Point 12: (1.0, -1.0, -1.0)\n');
fprintf('Same X,Y as P2, but Z=-1.0 (4 units above P2)\n');
[Stress, Strain] = TDstressHS(1.0, -1.0, -1.0, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('MATLAB Result: Exx = %.15e\n', Strain(1));
if isnan(Strain(1))
    fprintf('  -> NaN (singular)\n');
else
    fprintf('  -> Finite value\n');
end
fprintf('\n');

% Point 15: (1.0, -1.0, -8.0)
fprintf('------------------------------------------------------------\n');
fprintf('Point 15: (1.0, -1.0, -8.0)\n');
fprintf('Same X,Y as P2, but Z=-8.0 (3 units below P2)\n');
[Stress, Strain] = TDstressHS(1.0, -1.0, -8.0, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('MATLAB Result: Exx = %.15e\n', Strain(1));
if isnan(Strain(1))
    fprintf('  -> NaN (singular)\n');
else
    fprintf('  -> Finite value\n');
end
fprintf('\n');

fprintf('============================================================\n');
fprintf('INTERPRETATION:\n');
fprintf('These points lie on lines perpendicular to the triangle\n');
fprintf('passing through vertices (or close to them).\n');
fprintf('If MATLAB returns NaN: These are genuinely singular points\n');
fprintf('If MATLAB returns finite: Our Fortran edge detection is wrong\n');
fprintf('============================================================\n');
