% test_casez_matlab.m
% Test the casez_log case with specific points

clear all;
close all;

% Add path if needed (assuming TDstressHS.m is in current directory)

% Triangle vertices
p1 = [-1.0, -1.0, -5.0];
p2 = [1.0, -1.0, -5.0];
p3 = [-1.0, 1.0, -4.0];

% Test points
x = [-1.0/3.0, -1.0/3.0, -1.0/3.0, 7.0, -7.0, -1.0, -1.0, ...
     3.0, -3.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0]';
y = [-1.0/3.0, -1.0/3.0, -1.0/3.0, -1.0, -1.0, -3.0, 3.0, ...
     -3.0, 3.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0]';
z = [-3.0, -14.0/3.0, -6.0, -5.0, -5.0, -6.0, -3.0, ...
     -6.0, -3.0, -1.0, -1.0, -1.0, -8.0, -8.0, -8.0]';

% Slip components
ss = 1.0;   % Strike-slip
ds = -1.0;  % Dip-slip
ts = 2.0;   % Tensile-slip

% Elastic parameters
mu = 3.0e10;
lambda = 3.0e10;

% Calculate strains
n_points = length(x);
exx_results = zeros(n_points, 1);
trimode_results = zeros(n_points, 1);

fprintf('==============================================\n');
fprintf('Testing casez_log implementation (MATLAB)\n');
fprintf('==============================================\n');
fprintf('Triangle vertices:\n');
fprintf('  p1 = [%.4f, %.4f, %.4f]\n', p1);
fprintf('  p2 = [%.4f, %.4f, %.4f]\n', p2);
fprintf('  p3 = [%.4f, %.4f, %.4f]\n', p3);
fprintf('\n');
fprintf('Slip components: ss=%.1f, ds=%.1f, ts=%.1f\n', ss, ds, ts);
fprintf('Elastic params: mu=%.2e, lambda=%.2e\n', mu, lambda);
fprintf('\n');
fprintf('==============================================\n');
fprintf('Point#      x          y          z      trimode    e_xx\n');
fprintf('==============================================\n');

for i = 1:n_points
    % Call TDstressHS
    [stress, strain, ~] = TDstressHS(x(i), y(i), z(i), p1, p2, p3, ss, ds, ts, mu, lambda);

    exx_results(i) = strain(1);

    % Get trimode for this point
    trimode = trimodefinder(y(i), z(i), x(i), p1(2:3), p2(2:3), p3(2:3));
    trimode_results(i) = trimode;

    % Print results
    if isnan(strain(1))
        fprintf('%4d %10.4f %10.4f %10.4f %6d        NaN\n', ...
                i, x(i), y(i), z(i), trimode);
    else
        fprintf('%4d %10.4f %10.4f %10.4f %6d %15.6e\n', ...
                i, x(i), y(i), z(i), trimode, strain(1));
    end
end

fprintf('==============================================\n');
fprintf('\nSummary of trimode values:\n');
fprintf('  Points with trimode = +1 (casep): %d\n', sum(trimode_results == 1));
fprintf('  Points with trimode = -1 (casen): %d\n', sum(trimode_results == -1));
fprintf('  Points with trimode =  0 (casez): %d\n', sum(trimode_results == 0));
fprintf('==============================================\n');

% Identify which points have trimode == 0
casez_indices = find(trimode_results == 0);
if ~isempty(casez_indices)
    fprintf('\nPoints with casez_log (trimode=0):\n');
    for idx = casez_indices'
        fprintf('  Point %2d: x=%10.4f, y=%10.4f, z=%10.4f\n', ...
                idx, x(idx), y(idx), z(idx));
    end
end
