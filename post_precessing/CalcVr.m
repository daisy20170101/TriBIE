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
% calculate rupture speed of each SSE based on given slip rate threshold
% input:  slip rate output of TriBIE, 
% output: (along-strike) rupture speed 

%% load binary file of slip rate and T file
num_cell = 170496; % numnber of element in the mesh

folder = '/Volumes/seismology-1/Mexico/Mesh_all_depth/h90_N25_T2_May/'; % directory of data
twin = load(strcat([folder,'t_sse.txt']));
tfile = load(strcat([folder,'t_sse-h90_N25_T2.dat']));
tsse = tfile(:)*365;
dt = tsse(2:end) - tsse(1:end-1);

appendix = '-h95_N25_T1.dat';

jeve = 13; % # of episode

fin = strcat([folder,'data/dataSR',num2str(jeve),'.bin']); % load episode-based slip rate  file; binary file
fname = fopen(fin,'rb');
sr = fread(fname,[twin(jeve,3)*num_cell 1],'double');

tfile = load(strcat([folder,'t_sse',appendix]));
tsse = tfile(:)*365;

tt =  tsse(twin(jeve,1):twin(jeve,2))-tsse(twin(jeve,1));
dt =  tsse(twin(jeve,1)+1:twin(jeve,2))-tsse(twin(jeve,1):twin(jeve,2)-1);

fig_name = strcat([folder,'vr_',num2str(jeve),'.jpg']);

%% calculate rupture time at each element on the fault
Trup = zeros(num_cell,1)+1e9; % initial arrival time
mon = zeros(num_cell,1);
mu = 32e9; % shear modulus
yrs = 365*24*3600; % yrs to seconds
vc = log10(61*10/yrs/1d3); % threshold = log10(61*10/yrs/1d3)

for i = 1:twin(jeve,3)-1
    dd = [Trup,sr(i*num_cell-num_cell+1:num_cell*i)];
    
    dnum = find(dd(:,2) > vc & dd(:,1) > 1e8); 
    Trup(dnum) = tt(i);
    
    vcos = 0.5*(sr(i*num_cell-num_cell+1:i*num_cell)+sr(i*num_cell+1:i*num_cell+num_cell));
    dnum = find( vcos > vc); 
    mon(dnum) = mon(dnum) + 10.^vcos(dnum)*dt(i)*24*3600.*cellarea(dnum)*1e6*mu;
end

tnum = find(Trup<1e8); % find arrival time
ttotal = max(Trup(tnum)); % duration: maximum arrival time

mtotal = sum(mon);
m0 = 2/3*log10(mtotal*1e7)-10.7;

%% find element at 25 km depth
% ploting rupture propagation speed along the strike
% nnum = find(zc/1e3 > 25-0.02 & zc/1e3 < 25+0.02);
% [y_sort,ynum] = sort(yc(nnum));
nnew = load('depth25.txt');
y_sort = yc(nnew);

% calculte (along-strike) Vr based on gradient along 25 km depth-contour
trup_sort = Trup(nnew);
vr = (y_sort(2:end)-y_sort(1:end-1))/1e3./(trup_sort(2:end)-trup_sort(1:end-1));
vnum = find(abs(vr)<1e3);
vrtotal = mean(abs(vr(vnum))); % average rupture speed for each SSE

figure; 
set(gcf,'position',[100 100 750 300]);

subplot(1,2,2);
hold on; box on;
scatter(trup_sort(2:end),y_sort(2:end)/1e3,30,vr,'filled');
set(gca,'xlim',[0,500],'ylim',[-250,200]);
colormap(turbo);

h = colorbar;
caxis([0,4]);
h.Label.String ='Vr (km/day)';

xlabel('time (day)');
ylabel('y (km)');

subplot(1,2,1);
box on; hold on;
trisurf(trisse,Trup,'edgecolor','none');

h = colorbar;
caxis([0,500]); 
h.Label.String ='rupture time (day)';

saveas(gcf,fig_name,'jpeg');
