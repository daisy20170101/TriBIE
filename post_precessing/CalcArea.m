%%
% @file
% This file is part of TriBIE.
%
% @author Duo Li (https://www.geophysik.uni-muenchen.de/Members/dli)
%
% @section LICENSE
% Copyright (c) 2021, TriBIE Group
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% @section DESCRIPTION
% calculate area of each SSE based on given slip rate threshold
% input:  mesh file (vertice and topology)
% output: element area
%%
dir = '~/Documents/Mexico/Tri3DSSE/Mesh_all_depth/';
fin = fopen(strcat([dir,'400km_1km_smooth.gts']),'r'); % read in mesh file

nnum = textscan(fin,'%d %d %d\n',1);
nvex = nnum{1}; nedge = nnum{2}; ncell = nnum{3}; 
vex = textscan(fin,'%f %f %f\n',nvex);
cell = textscan(fin,'%d %d %d\n',ncell);
x = vex{1} ; y = vex{2}; z=vex{3};
x = x/1000; y = y/1000; z = z/1000;
n1 = cell{1}; n2 = cell{2}; n3 = cell{3};

x1 = x(n1);  y1=y(n1); z1=z(n1);
x2 = x(n2);  y2=y(n2); z2=z(n2);
x3 = x(n3);  y3=y(n3); z3=z(n3);

%% get nv(normal vector)nd ns (strike-vector) and element area.
nv = cross([x1-x2,y1-y2,z1-z2],[x1-x3,y1-y3,z1-z3]);
nv = nv./sqrt(nv(:,1).^2+nv(:,2).^2+nv(:,3).^2);

a = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2);
b = sqrt((x1-x3).^2+(y1-y3).^2+(z1-z3).^2);
c = sqrt((x3-x2).^2+(y3-y2).^2+(z3-z2).^2);
p = (a+b+c)/2;

cellarea = sqrt(p.*(p-a).*(p-b).*(p-c)); % area for triangular ele.

% centroid location
X = (x1+x2+x3)/3; 
Y = (y1+y2+y3)/3; 
Z = (z1+z2+z3)/3;

ns = zeros(length(nv(:,1)),3);
ns(:,1) = -sin(atan(nv(:,2)./nv(:,1)));
ns(:,2) = cos(atan(nv(:,2)./nv(:,1)));
nd = cross(nv,ns);

%% plot and save element area
figure; 
hold on;
box on;
trisurf([n1,n2,n3],x,y,z,cellarea,'edgecolor','none');

load('turbo');
colormap(turbo);
h = colorbar;
h.Label.String='element area (km^2)';
save('cellarea.txt', '-ascii', 'cellarea');
