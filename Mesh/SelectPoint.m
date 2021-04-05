%% written by D. Li

fin = fopen('400km_1km_smooth.gts','r');
% fin = fopen('triangular_mesh.gts','r'); % read in mesh file

nnum = textscan(fin,'%d %d %d\n',1);
nvex = nnum{1}; nedge = nnum{2}; ncell = nnum{3}; 
vex = textscan(fin,'%f %f %f\n',nvex);
cell = textscan(fin,'%d %d %d\n',ncell);
x = vex{1} ; 
y = vex{2}; 
z = vex{3};
x = x/1000; 
y = y/1000; 
z = z/1000;
n1 = cell{1}; 
n2 = cell{2}; 
n3 = cell{3};

x1 = x(n1);  y1=y(n1); z1=z(n1);
x2 = x(n2);  y2=y(n2); z2=z(n2);
x3 = x(n3);  y3=y(n3); z3=z(n3);

pp = [x1,y1,-z1];
% find points at depths
n1 = find(pp(:,3)> 9 & pp(:,3)< 10 & pp(:,2)>-132.0 & pp(:,2)<-130) ;
n2 = find(pp(:,3)> 19 & pp(:,3)< 20 & pp(:,2)>-132.0 & pp(:,2)<-130) ;
n3 = find(pp(:,3)> 29 & pp(:,3)< 30 & pp(:,2)>-132.0 & pp(:,2)<-130) ;
n4 = find(pp(:,3)> 39 & pp(:,3)< 40 & pp(:,2)>-132.0 & pp(:,2)<-130) ;

n5 = find(pp(:,3)> 19 & pp(:,3)< 20 & pp(:,2)>30.0 & pp(:,2)<32) ;
n6 = find(pp(:,3)> 19 & pp(:,3)< 20 & pp(:,2)>80.0 & pp(:,2)<82) ;
n7 = find(pp(:,3)> 19 & pp(:,3)< 20 & pp(:,2)>-82.0 & pp(:,2)<-80) ;
n8 = find(pp(:,3)> 19 & pp(:,3)< 20 & pp(:,2)>-32.0 & pp(:,2)<-30) ;

data1 = [n1(1),n2(1),n3(1),n4(1),n5(1),n6(1),n7(1),n8(1)];

save('ObvPoints.txt','-ascii','-int','data1');
%%
