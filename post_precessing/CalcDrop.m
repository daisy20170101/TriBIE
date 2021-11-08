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
% calculate stress drop of each SSE based on given slip rate threshold
% input:  shear traction output
% output: stress drop

%% load binary file
% num_cell = 170496;
% folder = 'h90_N25_T2_May/';
% twin = load(strcat([folder,'t_sse.txt']));
% tfile = load(strcat([folder,'t_sse-h90_N25_T2.dat']));
% tsse = tfile(:)*365;
% dt = tsse(2:end) - tsse(1:end-1);

% jeve = 13;

fin = strcat([folder,'data/dataTau',num2str(jeve),'.bin']); % load shear traction data
fname = fopen(fin,'rb');
tau = fread(fname,[twin(jeve,3)*num_cell 1],'double');

tt =  tsse(twin(jeve,1):twin(jeve,2))-tsse(twin(jeve,1));
figname = strcat([folder,'StressDrop_',num2str(jeve),'.jpg']);
%% calculate rupture time on the fault
Tpeak = zeros(num_cell,1);
Tresd = zeros(num_cell,1); 
Vpeak = zeros(num_cell,1); % peak slip rate

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

% dataout = [vptotal, vrtotal,tautotal, mtotal,m0,ttotal,jeve];
% save('scaling2_all.txt','dataout','-append','-ascii');

%%
figure; 
set(gcf,'position',[100 100 750 250]);
subplot(1,2,1);
hold on;box on;
trisurf(trisse,Tpeak,'edgecolor','none');
colormap(batlow);
h = colorbar;
h.Label.String ='shear traction (MPa)';
xlabel('longitude');
ylabel('latitude');
caxis([1.5 1.7]);

subplot(1,2,2);
hold on;box on;
trisurf(trisse,Tpeak-Tresd,'edgecolor','none');
colormap(vik);
h = colorbar;
caxis([0,0.3]); 
h.Label.String ='stress drop (MPa)';
xlabel('longitude');
ylabel('latitude');

saveas(gcf,figname,'jpeg');
