%% write geometry, unit is meter.
% generate triangular mesh of diameters approximating to  2000 m
fout = fopen('geometry_400slab_4.jou','w+');
fprintf(fout,'%s\n','${Units(si)}');
fprintf(fout,'%s\n','reset');
% fprintf(fout,'%s\n','undo off');

depth = [5 20 30 40 50 60];

for dd = 1:length(depth)
   fname = strcat('../../../2014SSE/CubitInput/csmooth2_',num2str(depth(dd)),'.txt');
   subgrd = load(fname);
   subgrd = subgrd(1:3:end,:);
   dnum = find(subgrd(:,2) < 170 & subgrd(:,2)> -230);
   subgrd = subgrd(dnum,:)*1e3;

%    subgrd(:,3) = -subgrd(:,3)-60e3 ;
   dep = (depth(dd)-4.9)*1e3; % depth should in positive.
   
   fprintf(fout,'%s%f%s%f%s%f%s\n','create vertex x {',subgrd(1,1),'} y {',subgrd(1,2),'} z {',dep,'}');  
   fprintf(fout,'%s%d%s\n','${idPtTopeS',dd,'=Id("vertex")}');
   
   for i = 2:length(subgrd(:,1))
    fprintf(fout,'%s%f%s%f%s%f%s\n','create vertex x {',subgrd(i,1),'} y {',subgrd(i,2),'} z {',dep,'}');  
   end
   fprintf(fout,'%s%d%s\n','${idPtTopeN',dd,'=Id("vertex")}');
%    fprintf(fout,'%s\n','${idPtTopeN=Id("vertex")}');
   
   fprintf(fout,'%s%d%s%d%s\n','create curve spline vertex {idPtTopeS',dd,'} to {idPtTopeN',dd,'}');
   fprintf(fout,'%s%d%s\n','curve {Id("curve")} name "curve',dd,'"');
   
%    if dd >1
%    fprintf(fout,'%s%d%s%d%s\n','create curve spline vertex {idPtTopeS',dd,'} {idPtTopeS',dd-1,'}'); 
%    fprintf(fout,'%s%d%s\n','curve {Id("curve")} name "edgeS',dd,'"');
%    fprintf(fout,'%s%d%s%d%s\n','create curve spline vertex {idPtTopeN',dd,'} {idPtTopeN',dd-1,'}');   
%    fprintf(fout,'%s%d%s\n','curve {Id("curve")} name "edgeN',dd,'"');
%    fprintf(fout,'%s%d%s%d%s%d%s%d\n','create surface curve curve',dd,' edgeN',dd,' curve',dd-1,' edgeS',dd);
%    fprintf(fout,'%s%d%s\n','surface {Id("surface")} name "surf',dd,'"');
%    end
end
% fprintf(fout,'%s%d %d %d %d %d %d %d %d %d %d %d\n','create surface curve ',1,2,3,4,5,6,7,8,9,10,11);
% fprintf(fout,'%s\n','create curve spline vertex {idPtTopeN1} {idPtTopeN2} {idPtTopeN3} {idPtTopeN4} {idPtTopeN5} {idPtTopeN6} {idPtTopeN7} {idPtTopeN8} {idPtTopeN9} {idPtTopeN10}');
% fprintf(fout,'%s\n','curve {Id("curve")} name "edgeN"');
% fprintf(fout,'%s\n','create curve spline vertex {idPtTopeS1} {idPtTopeS2} {idPtTopeS3} {idPtTopeS4} {idPtTopeS5} {idPtTopeS6} {idPtTopeS7} {idPtTopeS8} {idPtTopeS9} {idPtTopeS10}');
% fprintf(fout,'%s\n','curve {Id("curve")} name "edgeS"');
% fprintf(fout,'%s%d%s\n','create curve spline vertex {idPtTopeN1} {idPtTopeN',dd,'}');
% fprintf(fout,'%s\n','curve {Id("curve")} name "edgeN"');
% fprintf(fout,'%s%d%s\n','create curve spline vertex {idPtTopeS',dd,'} {idPtTopeS1}');
% fprintf(fout,'%s\n','curve {Id("curve")} name "edgeS"');
% fprintf(fout,'%s%d%s\n','create surface curve curve1 curve',dd,' edgeN edgeS');
% fprintf(fout,'%s\n','surface {Id("surface")} name "fault"');

% fprintf(fout,'%s\n','create surface vertex 100 120 130 on surface fault');
% fprintf(fout,'%s\n','delete vertex all');
% fprintf(fout,'%s\n','delete curve all');

%% create volume 
% fprintf(fout,'%s\n','create vertex x -300 y -800 z 0'); 
% fprintf(fout,'%s\n','${box1=Id("vertex")}');
% fprintf(fout,'%s\n','create vertex x 500 y -800 z 0'); 
% fprintf(fout,'%s\n','${box2=Id("vertex")}');
% fprintf(fout,'%s\n','create vertex x 500 y 300 z 0');  
% fprintf(fout,'%s\n','${box3=Id("vertex")}');
% fprintf(fout,'%s\n','create vertex x -300 y 300 z 0');  
% fprintf(fout,'%s\n','${box4=Id("vertex")}');
% fprintf(fout,'%s\n','create surface vertex {box1} {box2} {box3} {box4}'); 
% fprintf(fout,'%s\n','surface {Id("surface")} name "freesurf"');
% 
% fprintf(fout,'%s\n','sweep surface freesurf vector 0 0 -1 distance 200'); 
% fprintf(fout,'%s\n','volume {Id("volume")} name "faultbody"');
% fprintf(fout,'%s\n', 'Webcut volume faultbody sweep surface fault vector 0 0 1 Through_all');

% fprintf(fout,'%s\n','create vertex x -100 y -250 z 60');  
% fprintf(fout,'%s\n','create vertex x 300 y -250 z 60');  
% fprintf(fout,'%s\n','create vertex x -100 y 200 z 60');  
% fprintf(fout,'%s\n','create vertex x 300 y 200 z 60'); 
% fprintf(fout,'%s\n','delete surface fault');

fprintf(fout,'%s\n','create surface skin curve all');
fprintf(fout,'%s\n','delete curve all');
fprintf(fout,'%s\n','delete vertex all');
% fprintf(fout,'%s\n','export acis "surf_top.sat" overwrite');

% fprintf(fout,'%s\n','imprint all');
% fprintf(fout,'%s\n','merge all');
% fprintf(fout,'%s\n','rotate 45 about x');
% fprintf(fout,'%s\n','rotate 30 about z');
% fprintf(fout,'%s\n','save as "cascadia_mega.cub" overwrite');


% fprintf(fout,'%s\n','${dx=0.5}');
% fprintf(fout,'%s\n','surface fault scheme trimesh');
% fprintf(fout,'%s\n','mesh surface fault');


