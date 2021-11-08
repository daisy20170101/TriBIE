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
% extract source characteristics, Vr, Mo, stress drop, moment magnitude,  of each SSE based on given slip rate threshold
% input:  slip rate and shear traction 
% output: source charateristics

%% load binary file
num_cell = 170496;
folder = '/Volumes/seismology-1/Mexico/Mesh_all_depth/h90_N25_T2_May/';
twin = load(strcat([folder,'t_sse.txt']));
tfile = load(strcat([folder,'t_sse-h90_N25_T2.dat']));
tsse = tfile(:)*365;
dt = tsse(2:end) - tsse(1:end-1);

vc = log10(61*5/yrs/1d3); % threshold = log10(61*10/yrs/1d3); 10*vpl
mu = 32e9;

outname = strcat([folder,'source10vpl',appendix]);

   
for jeve = 21:24
    jeve
    
    Trup = zeros(num_cell,1)+1e9;
    mon = zeros(num_cell,1);
    Tpeak = zeros(num_cell,1);
    Tresd = zeros(num_cell,1); 
    Vpeak = zeros(num_cell,1); % peak slip rate
    
    fin = strcat([folder,'data/dataSR',num2str(jeve),'.bin']);
    fname = fopen(fin,'rb');
    sr = fread(fname,[twin(jeve,3)*num_cell 1],'double');

    tt =  tsse(twin(jeve,1):twin(jeve,2))-tsse(twin(jeve,1));
    dt =  tsse(twin(jeve,1)+1:twin(jeve,2))-tsse(twin(jeve,1):twin(jeve,2)-1);

    for i = 1:twin(jeve,3)-1
        dd = [Trup,sr(i*num_cell-num_cell+1:num_cell*i)];
        dnum = find(dd(:,2) > vc & dd(:,1) > 1e8); 
        Trup(dnum) = tt(i);

        vcos = 0.5*(sr(i*num_cell-num_cell+1:i*num_cell)+sr(i*num_cell+1:i*num_cell+num_cell));
        dnum = find( vcos > vc); 
        mon(dnum) = mon(dnum) + 10.^vcos(dnum)*dt(i)*24*3600.*cellarea(dnum)*1e6*mu;
    end
    tnum = find(Trup<1e8);
    ttotal = max(Trup(tnum));

    mtotal = sum(mon);
    m0 = 2/3*log10(mtotal*1e7)-10.7;

    % find element at 25 km depth
    % ploting rupture propagation speed along the strike
%     nnum = find(zc/1e3 > 25-0.02 & zc/1e3 < 25+0.02);
%     [y_sort,ynum] = sort(yc(nnum));

    trup_sort = Trup(nnew);
    vr = (y_sort(2:end)-y_sort(1:end-1))/1e3./(trup_sort(2:end)-trup_sort(1:end-1));% km/day
    vnum = find(abs(vr)<1e3);
    vrtotal = mean(abs(vr(vnum)));

    fin = strcat([folder,'data/dataTau',num2str(jeve),'.bin']);
    fname = fopen(fin,'rb');
    tau = fread(fname,[twin(jeve,3)*num_cell 1],'double');

    % calculate rupture time on the fault
   
    for i = 1:num_cell
       dtau = tau(i:num_cell:num_cell*twin(jeve,3));
       Tpeak(i) = max(dtau)/10;
       Tresd(i) = min(dtau)/10;

       dv = sr(i:num_cell:num_cell*twin(jeve,3));
       Vpeak(i) = max(dv);
    end

    vnum = find(Vpeak>vc);
    vptotal = max(Vpeak);
    tautotal = mean(Tpeak(vnum)-Tresd(vnum));

    dataout = [vptotal, vrtotal,tautotal, mtotal,m0,ttotal,jeve];
    save(outname,'dataout','-append','-ascii');
end