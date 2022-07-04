% @file
% This file is part of TriBIE.
%
% @author Duo Li dli@geophysik.uni-muenchen.de 
%
% @section LICENSE
% Copyright (c) 2022, TriBIE Group
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

%% read Inp file for mexico slab geometry. So Slow!! Use fortran ReadInp instead!!

% f1 = strcat('./mexicoslab_fine.inp');
% n_node = 654; % Test for Coulomb stress.
% n_ele = 1179;

% f1 = strcat('./400km_1km_smooth.inp');
fnod = strcat('./400_1km_node');
fele = strcat('./400_1km_ele');

n_node = 85884;
n_ele = 170534;

[ind1,xx,yy,zz] = textread(fnod,'%d,%f,%f,%f\n',n_node);
[ind2,n1,n2,n3] = textread(fele,'%d,%d,%d,%d\n',n_ele);

% figure; 
% hold on; box on;
for i = 1:n_ele
    ele1 = find(ind1 == n1(i));
    ele2 = find(ind1 == n2(i));
    ele3 = find(ind1 == n3(i));
    xl(i,:) = [xx(ele1),xx(ele2),xx(ele3),xx(ele1)];
    yl(i,:) = [yy(ele1),yy(ele2),yy(ele3),yy(ele1)];
    zl(i,:) = [zz(ele1),zz(ele2),zz(ele3),zz(ele1)];
    nn(i,1:3) = [ele1,ele2,ele3];
%     plot3(xl(i,:),yl(i,:),zl(i,:),'-k');
end

zz = -zz;

number = [n_node,0,n_ele];
data1 = [xx,yy,zz];

fout = fopen('fault_400km_1km_smooth.gts','w+'); %  dip angle may change

fprintf(fout,'%6d %6d %6d\n',number);

for i = 1:n_node
    fprintf(fout,'%f %f %f\n',data1(i,:));
    
end

for i = 1:n_ele
    fprintf(fout,'%6d %6d %6d\n',nn(i,:));
end
fclose(fout);
