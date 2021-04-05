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

data1 = [xx,yy,zz];
save 400km_1km_smooth.gts -ascii data1;
save 400km_1km_smooth.gts -ascii -append nn;