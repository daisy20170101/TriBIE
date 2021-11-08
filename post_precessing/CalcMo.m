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
% calculate moment magnitude of each SSE based on given slip rate threshold
% input:  slip rate data
% output: moment and equivalent moment magnitude

%% load binary file
% num_cell = 170496;
% folder = 'h90_N25_T2_May/';
% twin = load(strcat([folder,'t_sse.txt']));
% tfile = load(strcat([folder,'t_sse-h90_N25_T2.dat']));
% tsse = tfile(:)*365;

jeve = 13;

fin = strcat([folder,'data/dataSR',num2str(jeve),'.bin']);
fname = fopen(fin,'rb');
sr = fread(fname,[twin(jeve,3)*num_cell 1],'double');

tt =  tsse(twin(jeve,1):twin(jeve,2))-tsse(twin(jeve,1));
dt =  tsse(twin(jeve,1)+1:twin(jeve,2))-tsse(twin(jeve,1):twin(jeve,2)-1);

%% calculate rupture time on the fault
mon = zeros(num_cell,1);
mu = 32e9;
vc = log10(61*10/yrs/1d3); % threshold = log10(61*10/yrs/1d3)

for i = 1:twin(jeve,3)-1 % loop over duration of SSE
    vcos = 0.5*(sr(i*num_cell-num_cell+1:i*num_cell)+sr(i*num_cell+1:i*num_cell+num_cell));
    dnum = find( vcos > vc); 
    mon(dnum) = mon(dnum) + 10.^vcos(dnum)*dt(i)*24*3600.*cellarea(dnum)*1e6*mu;
end

mtotal= sum(mon);
m0 = 2/3*log10(mtotal*1e7)-10.7;

