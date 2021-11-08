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
% identify the onset and ending of each SSE based on given slip rate threshold

%% To identify event start and end time from t_sse/t-cos output file

folder = strcat('h98_N25_T1_2/');
tsse = load(strcat([folder, 't_sse-h95_N25_T1.dat']));
dtsse = tsse(2:end) - tsse(1:end-1);
threshold = 10/365; % identify time for SSE

T2 = find(dtsse > threshold); 
T1 = [1;T2(1:end-1)+1];

twin = [T1,T2,T2-T1+1];

outfile = strcat([folder,'/t_sse.txt']);

f1 = fopen(outfile,'w');

for i = 2:length(T1)
    fprintf(f1,'%i,%i,%i,%7.1f\n',T1(i),T2(i),T2(i)-T1(i)+1,tsse(T1(i)));
end