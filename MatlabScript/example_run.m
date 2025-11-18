% Example script for computing triangular Green's function matrix
%
% This script demonstrates how to use compute_trigreen_matrix.m
% to calculate the stiffness matrix for a triangular mesh.
%
% Author: TriBIE Development Team
% Date: 2025-11-18

%% Setup paths
% Add the NikkhooWalter2015 directory to MATLAB path
addpath('../NikkhooWalter2015');

% Check if TDstressHS.m exists
if ~exist('TDstressHS', 'file')
    error('TDstressHS.m not found. Please ensure ../NikkhooWalter2015 contains TDstressHS.m');
end

%% Example 1: Basic usage with default parameters
fprintf('\n=== Example 1: Basic usage ===\n');

% Define mesh file and number of processors
mesh_file = 'triangular_mesh.gts';  % Change to your mesh file
ncore = 4;                           % Number of processors

% Check if mesh file exists
if ~exist(mesh_file, 'file')
    fprintf('Warning: Mesh file %s not found.\n', mesh_file);
    fprintf('Please create a mesh file or update mesh_file variable.\n\n');
else
    % Compute with default parameters
    compute_trigreen_matrix(mesh_file, ncore);
end

%% Example 2: Custom material properties
fprintf('\n=== Example 2: Custom material properties ===\n');

if exist(mesh_file, 'file')
    % Specify custom material properties
    compute_trigreen_matrix(mesh_file, ncore, ...
        'mu', 32000, ...      % Shear modulus: 32 GPa
        'nu', 0.28, ...       % Poisson's ratio: 0.28
        'slip_ss', -1.0, ...  % Strike-slip: -1.0
        'slip_ds', 0.0, ...   % Dip-slip: 0.0
        'slip_ts', 0.0);      % Tensile-slip: 0.0
end

%% Example 3: Different slip components
fprintf('\n=== Example 3: Pure dip-slip ===\n');

if exist(mesh_file, 'file')
    % Pure dip-slip configuration
    compute_trigreen_matrix(mesh_file, ncore, ...
        'slip_ss', 0.0, ...   % No strike-slip
        'slip_ds', 1.0, ...   % Pure dip-slip
        'slip_ts', 0.0);      % No tensile-slip
end

%% Example 4: Specify output directory
fprintf('\n=== Example 4: Custom output directory ===\n');

if exist(mesh_file, 'file')
    output_dir = './output_trigreen';

    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Compute with custom output directory
    compute_trigreen_matrix(mesh_file, ncore, ...
        'output_dir', output_dir);
end

%% Example 5: Single processor
fprintf('\n=== Example 5: Single processor computation ===\n');

if exist(mesh_file, 'file')
    % Run on single processor
    compute_trigreen_matrix(mesh_file, 1);
end

%% Example 6: Parallel configuration (many processors)
fprintf('\n=== Example 6: Many processors ===\n');

if exist(mesh_file, 'file')
    % Distribute across 16 processors
    compute_trigreen_matrix(mesh_file, 16);
end

fprintf('\n=== Examples completed ===\n');
fprintf('Check output directory for trigreen_*.bin files\n\n');
