% Test Points 8 and 9 with detailed contribution breakdown
% Compare with Fortran output
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

fprintf('Testing Points 8 and 9 - Detailed Contributions\n');
fprintf('===========================================================\n\n');

% Point 8
X = 3.0;
Y = -3.0;
Z = -6.0;

fprintf('Point 8: (%.1f, %.1f, %.1f)\n', X, Y, Z);
fprintf('-----------------------------------------------------------\n');

% Main dislocation (full-space)
[StsMS, StrMS] = TDstressFS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Main Dislocation:\n');
fprintf('  Exx = %.15e\n', StrMS(1));
fprintf('  Eyy = %.15e\n', StrMS(2));
fprintf('  Ezz = %.15e\n', StrMS(3));

% Harmonic function
[StsFSC, StrFSC] = TDstress_HarFunc(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Harmonic Function:\n');
fprintf('  Exx = %.15e\n', StrFSC(1));
fprintf('  Eyy = %.15e\n', StrFSC(2));
fprintf('  Ezz = %.15e\n', StrFSC(3));

% Image dislocation (flip z-coordinates)
P1_img = P1; P1_img(3) = -P1_img(3);
P2_img = P2; P2_img(3) = -P2_img(3);
P3_img = P3; P3_img(3) = -P3_img(3);
[StsIS, StrIS] = TDstressFS(X, Y, Z, P1_img, P2_img, P3_img, Ss, Ds, Ts, mu, lambda);
fprintf('Image Dislocation:\n');
fprintf('  Exx = %.15e\n', StrIS(1));
fprintf('  Eyy = %.15e\n', StrIS(2));
fprintf('  Ezz = %.15e\n', StrIS(3));

% Total
Total_Exx = StrMS(1) + StrFSC(1) + StrIS(1);
Total_Eyy = StrMS(2) + StrFSC(2) + StrIS(2);
Total_Ezz = StrMS(3) + StrFSC(3) + StrIS(3);
fprintf('Total (Main + Harmonic + Image):\n');
fprintf('  Exx = %.15e\n', Total_Exx);
fprintf('  Eyy = %.15e\n', Total_Eyy);
fprintf('  Ezz = %.15e\n', Total_Ezz);
fprintf('\n');

% Point 9
X = -3.0;
Y = 3.0;
Z = -3.0;

fprintf('Point 9: (%.1f, %.1f, %.1f)\n', X, Y, Z);
fprintf('-----------------------------------------------------------\n');

% Main dislocation (full-space)
[StsMS, StrMS] = TDstressFS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Main Dislocation:\n');
fprintf('  Exx = %.15e\n', StrMS(1));
fprintf('  Eyy = %.15e\n', StrMS(2));
fprintf('  Ezz = %.15e\n', StrMS(3));

% Harmonic function
[StsFSC, StrFSC] = TDstress_HarFunc(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lambda);
fprintf('Harmonic Function:\n');
fprintf('  Exx = %.15e\n', StrFSC(1));
fprintf('  Eyy = %.15e\n', StrFSC(2));
fprintf('  Ezz = %.15e\n', StrFSC(3));

% Image dislocation (flip z-coordinates)
P1_img = P1; P1_img(3) = -P1_img(3);
P2_img = P2; P2_img(3) = -P2_img(3);
P3_img = P3; P3_img(3) = -P3_img(3);
[StsIS, StrIS] = TDstressFS(X, Y, Z, P1_img, P2_img, P3_img, Ss, Ds, Ts, mu, lambda);
fprintf('Image Dislocation:\n');
fprintf('  Exx = %.15e\n', StrIS(1));
fprintf('  Eyy = %.15e\n', StrIS(2));
fprintf('  Ezz = %.15e\n', StrIS(3));

% Total
Total_Exx = StrMS(1) + StrFSC(1) + StrIS(1);
Total_Eyy = StrMS(2) + StrFSC(2) + StrIS(2);
Total_Ezz = StrMS(3) + StrFSC(3) + StrIS(3);
fprintf('Total (Main + Harmonic + Image):\n');
fprintf('  Exx = %.15e\n', Total_Exx);
fprintf('  Eyy = %.15e\n', Total_Eyy);
fprintf('  Ezz = %.15e\n', Total_Ezz);
fprintf('\n');

fprintf('===========================================================\n');
fprintf('Expected Total Values:\n');
fprintf('  Point 8: Exx = 7.064e-4\n');
fprintf('  Point 9: Exx = 2.113e-4\n');
fprintf('===========================================================\n');
