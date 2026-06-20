% debug_matlab.m
% MATLAB debug script to generate reference values for comparison

% Test parameters (matching Fortran debug program)
X = -1/3;
Y = -1/3;
Z = -14/3;

P1 = [-1, -1, -5];
P2 = [1, -1, -5];
P3 = [-1, 1, -4];

Ss = 1.0;   % Strike-slip
Ds = -1.0;  % Dip-slip
Ts = 2.0;   % Tensile-slip

mu = 3.0e10;      % Shear modulus
lambda = 3.0e10;  % Lame's first parameter

fprintf('=== MATLAB DEBUG: Nikkhoo & Walter (2015) ===\n');
fprintf('Input Parameters:\n');
fprintf('  Calculation point: (%.6f, %.6f, %.6f)\n', X, Y, Z);
fprintf('  Triangle vertices:\n');
fprintf('    P1 = (%.6f, %.6f, %.6f)\n', P1(1), P1(2), P1(3));
fprintf('    P2 = (%.6f, %.6f, %.6f)\n', P2(1), P2(2), P2(3));
fprintf('    P3 = (%.6f, %.6f, %.6f)\n', P3(1), P3(2), P3(3));
fprintf('  Slip components: SS = %.6f, DS = %.6f, TS = %.6f\n', Ss, Ds, Ts);
fprintf('  Elastic parameters: mu = %.6e, lambda = %.6e\n', mu, lambda);

% Calculate Poisson's ratio
nu = 1/(1+lambda/mu)/2;
fprintf('  Poisson ratio: nu = %.6f\n', nu);

% Slip vector components
bx = Ts;  % Tensile-slip
by = Ss;  % Strike-slip
bz = Ds;  % Dip-slip
fprintf('  Slip vector: (%.6f, %.6f, %.6f)\n', bx, by, bz);
fprintf('\n');

% Calculate unit vectors
ey = [0, 1, 0];
ez = [0, 0, 1];

% Normal vector
vnorm = cross(P2 - P1, P3 - P1);
vnorm = vnorm / norm(vnorm);

% Strike vector
vstrike = cross(ez, vnorm);
if norm(vstrike) < eps
    vstrike = ey * vnorm(3);
    if P1(3) > 0
        vstrike = -vstrike;
    end
end
vstrike = vstrike / norm(vstrike);

% Dip vector
vdip = cross(vnorm, vstrike);

fprintf('=== Coordinate System Vectors ===\n');
fprintf('  vnorm = (%.6f, %.6f, %.6f)\n', vnorm(1), vnorm(2), vnorm(3));
fprintf('  vstrike = (%.6f, %.6f, %.6f)\n', vstrike(1), vstrike(2), vstrike(3));
fprintf('  vdip = (%.6f, %.6f, %.6f)\n', vdip(1), vdip(2), vdip(3));
fprintf('\n');

% Transformation matrix (transpose as in MATLAB)
A = [vnorm; vstrike; vdip];

fprintf('=== Transformation Matrix A (EFCS to TDCS) ===\n');
fprintf('  A(1,:) = (%.6f, %.6f, %.6f)\n', A(1,1), A(1,2), A(1,3));
fprintf('  A(2,:) = (%.6f, %.6f, %.6f)\n', A(2,1), A(2,2), A(2,3));
fprintf('  A(3,:) = (%.6f, %.6f, %.6f)\n', A(3,1), A(3,2), A(3,3));
fprintf('\n');

% Transform coordinates to TDCS
p1 = zeros(3,1);
p2 = zeros(3,1);
p3 = zeros(3,1);

[x,y,z] = CoordTrans(X-P2(1),Y-P2(2),Z-P2(3),A);
[p1(1),p1(2),p1(3)] = CoordTrans(P1(1)-P2(1),P1(2)-P2(2),P1(3)-P2(3),A);
[p3(1),p3(2),p3(3)] = CoordTrans(P3(1)-P2(1),P3(2)-P2(2),P3(3)-P2(3),A);

fprintf('=== TDCS Coordinates ===\n');
fprintf('  Calculation point: (%.6f, %.6f, %.6f)\n', x, y, z);
fprintf('  Triangle vertices:\n');
fprintf('    p1 = (%.6f, %.6f, %.6f)\n', p1(1), p1(2), p1(3));
fprintf('    p2 = (%.6f, %.6f, %.6f)\n', p2(1), p2(2), p2(3));
fprintf('    p3 = (%.6f, %.6f, %.6f)\n', p3(1), p3(2), p3(3));
fprintf('\n');

% Calculate unit vectors along TD sides
e12 = (p2-p1)/norm(p2-p1);
e13 = (p3-p1)/norm(p3-p1);
e23 = (p3-p2)/norm(p3-p2);

fprintf('=== Unit Vectors Along TD Sides ===\n');
fprintf('  e12 = (%.6f, %.6f, %.6f)\n', e12(1), e12(2), e12(3));
fprintf('  e13 = (%.6f, %.6f, %.6f)\n', e13(1), e13(2), e13(3));
fprintf('  e23 = (%.6f, %.6f, %.6f)\n', e23(1), e23(2), e23(3));
fprintf('\n');

% Calculate angles
A_angle = acos(e12'*e13);
B_angle = acos(-e12'*e23);
C_angle = acos(e23'*e13);

fprintf('=== Triangle Angles ===\n');
fprintf('  A_angle = %.6f rad = %.6f deg\n', A_angle, A_angle * 180 / pi);
fprintf('  B_angle = %.6f rad = %.6f deg\n', B_angle, B_angle * 180 / pi);
fprintf('  C_angle = %.6f rad = %.6f deg\n', C_angle, C_angle * 180 / pi);
fprintf('\n');

% Determine configuration
Trimode = trimodefinder(y,z,x,p1(2:3),p2(2:3),p3(2:3));

fprintf('=== Configuration ===\n');
fprintf('  trimode = %d\n', Trimode);
fprintf('  casep_log = %d\n', Trimode==1);
fprintf('  casen_log = %d\n', Trimode==-1);
fprintf('  casez_log = %d\n', Trimode==0);
fprintf('\n');

% Calculate individual contributions for debugging
fprintf('=== Calling TDstressHS with detailed breakdown ===\n');

% Calculate main dislocation contribution
[Stress_MS,Strain_MS] = TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);

fprintf('=== Main Dislocation Contribution ===\n');
fprintf('Stress: Sxx=%.6e Syy=%.6e Szz=%.6e Sxy=%.6e Sxz=%.6e Syz=%.6e\n', ...
        Stress_MS(1), Stress_MS(2), Stress_MS(3), Stress_MS(4), Stress_MS(5), Stress_MS(6));
fprintf('Strain: Exx=%.6e Eyy=%.6e Ezz=%.6e Exy=%.6e Exz=%.6e Eyz=%.6e\n', ...
        Strain_MS(1), Strain_MS(2), Strain_MS(3), Strain_MS(4), Strain_MS(5), Strain_MS(6));
fprintf('\n');

% Calculate harmonic function contribution
[Stress_FSC,Strain_FSC] = TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);

fprintf('=== Harmonic Function Contribution ===\n');
fprintf('Stress: Sxx=%.6e Syy=%.6e Szz=%.6e Sxy=%.6e Sxz=%.6e Syz=%.6e\n', ...
        Stress_FSC(1), Stress_FSC(2), Stress_FSC(3), Stress_FSC(4), Stress_FSC(5), Stress_FSC(6));
fprintf('Strain: Exx=%.6e Eyy=%.6e Ezz=%.6e Exy=%.6e Exz=%.6e Eyz=%.6e\n', ...
        Strain_FSC(1), Strain_FSC(2), Strain_FSC(3), Strain_FSC(4), Strain_FSC(5), Strain_FSC(6));
fprintf('\n');

% Calculate image dislocation contribution
P1_img = P1; P2_img = P2; P3_img = P3;
P1_img(3) = -P1_img(3);
P2_img(3) = -P2_img(3);
P3_img(3) = -P3_img(3);

[Stress_IS,Strain_IS] = TDstressFS(X,Y,Z,P1_img,P2_img,P3_img,Ss,Ds,Ts,mu,lambda);

% Special case for surface elements
if abs(P1_img(3)) < eps && abs(P2_img(3)) < eps && abs(P3_img(3)) < eps
    Stress_IS(5) = -Stress_IS(5);  % xz component
    Stress_IS(6) = -Stress_IS(6);  % yz component
    Strain_IS(5) = -Strain_IS(5);  % xz component
    Strain_IS(6) = -Strain_IS(6);  % yz component
    fprintf('Applied surface element correction\n');
end

fprintf('=== Image Dislocation Contribution ===\n');
fprintf('Stress: Sxx=%.6e Syy=%.6e Szz=%.6e Sxy=%.6e Sxz=%.6e Syz=%.6e\n', ...
        Stress_IS(1), Stress_IS(2), Stress_IS(3), Stress_IS(4), Stress_IS(5), Stress_IS(6));
fprintf('Strain: Exx=%.6e Eyy=%.6e Ezz=%.6e Exy=%.6e Exz=%.6e Eyz=%.6e\n', ...
        Strain_IS(1), Strain_IS(2), Strain_IS(3), Strain_IS(4), Strain_IS(5), Strain_IS(6));
fprintf('\n');

% Calculate total results
Stress_total = Stress_MS + Stress_IS + Stress_FSC;
Strain_total = Strain_MS + Strain_IS + Strain_FSC;

fprintf('=== Total Results ===\n');
fprintf('Stress: Sxx=%.6e Syy=%.6e Szz=%.6e Sxy=%.6e Sxz=%.6e Syz=%.6e\n', ...
        Stress_total(1), Stress_total(2), Stress_total(3), Stress_total(4), Stress_total(5), Stress_total(6));
fprintf('Strain: Exx=%.6e Eyy=%.6e Ezz=%.6e Exy=%.6e Exz=%.6e Eyz=%.6e\n', ...
        Strain_total(1), Strain_total(2), Strain_total(3), Strain_total(4), Strain_total(5), Strain_total(6));
fprintf('\n');

% Calculate final results using the original function
[Stress,Strain] = TDstressHS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);

fprintf('=== Final Results (from TDstressHS) ===\n');
fprintf('Stress tensor:\n');
fprintf('  Sxx = %.6e\n', Stress(1));
fprintf('  Syy = %.6e\n', Stress(2));
fprintf('  Szz = %.6e\n', Stress(3));
fprintf('  Sxy = %.6e\n', Stress(4));
fprintf('  Sxz = %.6e\n', Stress(5));
fprintf('  Syz = %.6e\n', Stress(6));
fprintf('\n');
fprintf('Strain tensor:\n');
fprintf('  Exx = %.6e\n', Strain(1));
fprintf('  Eyy = %.6e\n', Strain(2));
fprintf('  Ezz = %.6e\n', Strain(3));
fprintf('  Exy = %.6e\n', Strain(4));
fprintf('  Exz = %.6e\n', Strain(5));
fprintf('  Eyz = %.6e\n', Strain(6));
fprintf('\n');

function [X1,X2,X3]=CoordTrans(x1,x2,x3,A)
% CoordTrans transforms the coordinates of the vectors, from
% x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
% transformation matrix, whose columns e1,e2 and e3 are the unit base 
% vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
% in X1X2X3. The transpose of A (i.e., A') will transform the coordinates 
% from X1X2X3 into x1x2x3.

x1 = x1(:);
x2 = x2(:);
x3 = x3(:);

% Debug: Check dimensions
fprintf('Debug CoordTrans: A size = [%d, %d], x1 size = [%d, %d], x2 size = [%d, %d], x3 size = [%d, %d]\n', ...
        size(A,1), size(A,2), size(x1,1), size(x1,2), size(x2,1), size(x2,2), size(x3,1), size(x3,2));

% Handle both scalar and vector inputs
if length(x1) == 1
    % Scalar case - A should be 3x3, [x1;x2;x3] is 3x1, result is 3x1
    input_vec = [x1; x2; x3];
    fprintf('Debug: input_vec size = [%d, %d]\n', size(input_vec,1), size(input_vec,2));
    r = A * input_vec;
    X1 = r(1);
    X2 = r(2);
    X3 = r(3);
else
    % Vector case - A should be 3x3, [x1';x2';x3'] is 3xN, result is 3xN
    input_mat = [x1'; x2'; x3'];
    fprintf('Debug: input_mat size = [%d, %d]\n', size(input_mat,1), size(input_mat,2));
    r = A * input_mat;
    X1 = r(1, :)';
    X2 = r(2, :)';
    X3 = r(3, :)';
end

end

function [trimode]=trimodefinder(x,y,z,p1,p2,p3)
% trimodefinder calculates the normalized barycentric coordinates of 
% the points with respect to the TD vertices and specifies the appropriate
% artefact-free configuration of the angular dislocations for the 
% calculations. The input matrices x, y and z share the same size and
% correspond to the y, z and x coordinates in the TDCS, respectively. p1,
% p2 and p3 are two-component matrices representing the y and z coordinates
% of the TD vertices in the TDCS, respectively.
% The components of the output (trimode) corresponding to each calculation 
% points, are 1 for the first configuration, -1 for the second 
% configuration and 0 for the calculation point that lie on the TD sides.

x = x(:);
y = y(:);
z = z(:);

a = ((p2(2)-p3(2)).*(x-p3(1))+(p3(1)-p2(1)).*(y-p3(2)))./...
    ((p2(2)-p3(2)).*(p1(1)-p3(1))+(p3(1)-p2(1)).*(p1(2)-p3(2)));
b = ((p3(2)-p1(2)).*(x-p3(1))+(p1(1)-p3(1)).*(y-p3(2)))./...
    ((p2(2)-p3(2)).*(p1(1)-p3(1))+(p3(1)-p2(1)).*(p1(2)-p3(2)));
c = 1-a-b;

trimode = ones(length(x),1);
trimode(a<=0 & b>c & c>a) = -1;
trimode(b<=0 & c>a & a>b) = -1;
trimode(c<=0 & a>b & b>c) = -1;
trimode(a==0 & b>=0 & c>=0) = 0;
trimode(a>=0 & b==0 & c>=0) = 0;
trimode(a>=0 & b>=0 & c==0) = 0;
trimode(trimode==0 & z~=0) = 1;
end

function [Stress,Strain]=TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda)
% TDstress_HarFunc calculates the harmonic function contribution to the
% strains and stresses associated with a triangular dislocation in a 
% half-space. The function cancels the surface normal tractions induced by 
% the main and image dislocations.

bx = Ts; % Tensile-slip
by = Ss; % Strike-slip
bz = Ds; % Dip-slip

% Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
% as an exception, if the normal vector points upward, the strike and dip 
% vectors point Northward and Westward, whereas if the normal vector points
% downward, the strike and dip vectors point Southward and Westward, 
% respectively.
Vnorm = cross(P2-P1,P3-P1);
Vnorm = Vnorm(:);  % Ensure column vector
Vnorm = Vnorm/norm(Vnorm);

eY = [0 1 0]';
eZ = [0 0 1]';
Vstrike = cross(eZ,Vnorm);

if norm(Vstrike)==0
    Vstrike = eY*Vnorm(3);
end
Vstrike = Vstrike(:);  % Ensure column vector
Vstrike = Vstrike/norm(Vstrike);
Vdip = cross(Vnorm,Vstrike);
Vdip = Vdip(:);  % Ensure column vector

% Debug: Check dimensions before creating A
fprintf('Debug TDstress_HarFunc: Vnorm size = [%d, %d], Vstrike size = [%d, %d], Vdip size = [%d, %d]\n', ...
        size(Vnorm,1), size(Vnorm,2), size(Vstrike,1), size(Vstrike,2), size(Vdip,1), size(Vdip,2));

% Transform slip vector components from TDCS into EFCS
A = [Vnorm Vstrike Vdip];
fprintf('Debug TDstress_HarFunc: A size = [%d, %d]\n', size(A,1), size(A,2));
[bX,bY,bZ] = CoordTrans(bx,by,bz,A);

% Calculate contribution of angular dislocation pair on each TD side 
[Stress1,Strain1] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P1,P2,mu,lambda); % P1P2
[Stress2,Strain2] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P2,P3,mu,lambda); % P2P3
[Stress3,Strain3] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P3,P1,mu,lambda); % P3P1

% Calculate total harmonic function contribution to strains and stresses
Stress = Stress1+Stress2+Stress3;
Strain = Strain1+Strain2+Strain3;
end

function [Stress,Strain]=AngSetupFSC_S(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda)
% AngSetupFSC_S calculates the Free Surface Correction to strains and 
% stresses associated with angular dislocation pair on each TD side.

nu = 1/(1+lambda/mu)/2; % Poisson's ratio

% Calculate TD side vector and the angle of the angular dislocation pair
SideVec = PB-PA;
SideVec = SideVec(:);  % Ensure column vector
eZ = [0 0 1]';
beta = acos(-SideVec'*eZ/norm(SideVec));

if abs(beta)<eps || abs(pi-beta)<eps
    Stress = zeros(length(X),6);
    Strain = zeros(length(X),6);
else
    ey1 = [SideVec(1:2);0];
    ey1 = ey1/norm(ey1);
    ey3 = -eZ;
    ey2 = cross(ey3,ey1);
    A = [ey1,ey2,ey3]; % Transformation matrix
    
    % Transform coordinates from EFCS to the first ADCS
    [y1A,y2A,y3A] = CoordTrans(X-PA(1),Y-PA(2),Z-PA(3),A);
    % Transform coordinates from EFCS to the second ADCS
    [y1AB,y2AB,y3AB] = CoordTrans(SideVec(1),SideVec(2),SideVec(3),A);
    y1B = y1A-y1AB;
    y2B = y2A-y2AB;
    y3B = y3A-y3AB;
    
    % Transform slip vector components from EFCS to ADCS
    [b1,b2,b3] = CoordTrans(bX,bY,bZ,A);
    
    % Determine the best arteact-free configuration for the calculation
    % points near the free furface
    I = (beta*y1A)>=0;
    
    % For singularities at surface
    v11A = zeros(length(X),1);
    v22A = zeros(length(X),1);
    v33A = zeros(length(X),1);
    v12A = zeros(length(X),1);
    v13A = zeros(length(X),1);
    v23A = zeros(length(X),1);
    
    v11B = zeros(length(X),1);
    v22B = zeros(length(X),1);
    v33B = zeros(length(X),1);
    v12B = zeros(length(X),1);
    v13B = zeros(length(X),1);
    v23B = zeros(length(X),1);
    
    % Configuration I
    [v11A(I),v22A(I),v33A(I),v12A(I),v13A(I),v23A(I)] = ...
        AngDisStrainFSC(-y1A(I),-y2A(I),y3A(I),...
        pi-beta,-b1,-b2,b3,nu,-PA(3));
    v13A(I) = -v13A(I);
    v23A(I) = -v23A(I);
    
    [v11B(I),v22B(I),v33B(I),v12B(I),v13B(I),v23B(I)] = ...
        AngDisStrainFSC(-y1B(I),-y2B(I),y3B(I),...
        pi-beta,-b1,-b2,b3,nu,-PB(3));
    v13B(I) = -v13B(I);
    v23B(I) = -v23B(I);
    
    % Configuration II
    [v11A(~I),v22A(~I),v33A(~I),v12A(~I),v13A(~I),v23A(~I)] = ...
        AngDisStrainFSC(y1A(~I),y2A(~I),y3A(~I),...
        beta,b1,b2,b3,nu,-PA(3));
    
    [v11B(~I),v22B(~I),v33B(~I),v12B(~I),v13B(~I),v23B(~I)] = ...
        AngDisStrainFSC(y1B(~I),y2B(~I),y3B(~I),...
        beta,b1,b2,b3,nu,-PB(3));
    
    % Calculate total Free Surface Correction to strains in ADCS
    v11 = v11B-v11A;
    v22 = v22B-v22A;
    v33 = v33B-v33A;
    v12 = v12B-v12A;
    v13 = v13B-v13A;
    v23 = v23B-v23A;
    
    % Transform total Free Surface Correction to strains from ADCS to EFCS
    [Exx,Eyy,Ezz,Exy,Exz,Eyz] = TensTrans(v11,v22,v33,v12,v13,v23,A');
    
    % Calculate total Free Surface Correction to stresses in EFCS
    Sxx = 2*mu*Exx+lambda*(Exx+Eyy+Ezz);
    Syy = 2*mu*Eyy+lambda*(Exx+Eyy+Ezz);
    Szz = 2*mu*Ezz+lambda*(Exx+Eyy+Ezz);
    Sxy = 2*mu*Exy;
    Sxz = 2*mu*Exz;
    Syz = 2*mu*Eyz;
    
    Strain = [Exx,Eyy,Ezz,Exy,Exz,Eyz];
    Stress = [Sxx,Syy,Szz,Sxy,Sxz,Syz];
end

end

function [Exx,Eyy,Ezz,Exy,Exz,Eyz]=AngDisStrainFSC(x,y,z,alpha,bx,by,bz,nu,a)
% AngDisStrainFSC calculates the strains associated with an angular 
% dislocation in an elastic half-space with free surface correction.

% This is a simplified placeholder - you'll need to implement the full function
% based on the original MATLAB code
Exx = zeros(size(x));
Eyy = zeros(size(x));
Ezz = zeros(size(x));
Exy = zeros(size(x));
Exz = zeros(size(x));
Eyz = zeros(size(x));
end

function [Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2]=TensTrans(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,A)
% TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate
% system to x2y2z2 coordinate system. "A" is the transformation matrix, 
% whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The 
% coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
% of A (i.e., A') does the transformation from x2y2z2 into x1y1z1.

Txx2 = A(1)^2*Txx1+2*A(1)*A(4)*Txy1+2*A(1)*A(7)*Txz1+2*A(4)*A(7)*Tyz1+...
    A(4)^2*Tyy1+A(7)^2*Tzz1;

Tyy2 = A(2)^2*Txx1+2*A(2)*A(5)*Txy1+2*A(2)*A(8)*Txz1+2*A(5)*A(8)*Tyz1+...
    A(5)^2*Tyy1+A(8)^2*Tzz1;

Tzz2 = A(3)^2*Txx1+2*A(3)*A(6)*Txy1+2*A(3)*A(9)*Txz1+2*A(6)*A(9)*Tyz1+...
    A(6)^2*Tyy1+A(9)^2*Tzz1;

Txy2 = A(1)*A(2)*Txx1+(A(1)*A(5)+A(2)*A(4))*Txy1+(A(1)*A(8)+...
    A(2)*A(7))*Txz1+(A(8)*A(4)+A(7)*A(5))*Tyz1+A(5)*A(4)*Tyy1+...
    A(7)*A(8)*Tzz1;

Txz2 = A(1)*A(3)*Txx1+(A(1)*A(6)+A(3)*A(4))*Txy1+(A(1)*A(9)+...
    A(3)*A(7))*Txz1+(A(9)*A(4)+A(7)*A(6))*Tyz1+A(6)*A(4)*Tyy1+...
    A(7)*A(9)*Tzz1;

Tyz2 = A(2)*A(3)*Txx1+(A(3)*A(5)+A(2)*A(6))*Txy1+(A(3)*A(8)+...
    A(2)*A(9))*Txz1+(A(8)*A(6)+A(9)*A(5))*Tyz1+A(5)*A(6)*Tyy1+...
    A(8)*A(9)*Tzz1;
end