% read hdf5 

h5file = strcat('output/step07a-fault-slab.h5');

h5disp(h5file,'/vertex_fields/traction_change');
h5disp(h5file,'/topology/cells');
h5disp(h5file,'/geometry/vertices');

connect = h5read(h5file,'/topology/cells/');
connect = connect +1 ;
xyz = h5read(h5file,'/geometry/vertices');
xyz = xyz';
connect = connect';

tri = triangulation(connect, xyz(:,1),xyz(:,2),xyz(:,3));

%%
figure;
trisurf(tri);



