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
%% write geometry, unit is meter.
% generate triangular mesh of diameters approximating to  2000 m
fout = fopen('geometry_400slab_4.jou','w+');
fprintf(fout,'%s\n','${Units(si)}');
fprintf(fout,'%s\n','reset');
% fprintf(fout,'%s\n','undo off');

% setup contour files' folder, prefix and depth steps
prefix = strcat('../../../2014SSE/CubitInput/csmooth2_'); % folder + contour files' prefix
depth = [5 20 30 40 50 60];

for dd = 1:length(depth)
   fname = strcat(prefix,num2str(depth(dd)),'.txt');
   subgrd = load(fname);
   subgrd = subgrd(1:3:end,:);

   dnum = find(subgrd(:,2) < 170.0 & subgrd(:,2)> -230.0); % if use subdomain from original contour file

   subgrd = subgrd(dnum,:)*1e3;

   dep = (depth(dd))*1e3; % pay attention depth should in positive.
   
   fprintf(fout,'%s%f%s%f%s%f%s\n','create vertex x {',subgrd(1,1),'} y {',subgrd(1,2),'} z {',dep,'}');  
   fprintf(fout,'%s%d%s\n','${idPtTopeS',dd,'=Id("vertex")}');
   
   for i = 2:length(subgrd(:,1))
    fprintf(fout,'%s%f%s%f%s%f%s\n','create vertex x {',subgrd(i,1),'} y {',subgrd(i,2),'} z {',dep,'}');  
   end
   fprintf(fout,'%s%d%s\n','${idPtTopeN',dd,'=Id("vertex")}');
%    fprintf(fout,'%s\n','${idPtTopeN=Id("vertex")}');
   
   fprintf(fout,'%s%d%s%d%s\n','create curve spline vertex {idPtTopeS',dd,'} to {idPtTopeN',dd,'}');
   fprintf(fout,'%s%d%s\n','curve {Id("curve")} name "curve',dd,'"');
   
end

%% create volume 

fprintf(fout,'%s\n','create surface skin curve all');
fprintf(fout,'%s\n','delete curve all');
fprintf(fout,'%s\n','delete vertex all');


