function vis20Examples
tic
figGenSFC3DCase('./data/testVols/nucleon32.nrrd');
toc

tic
figGenSFC2DCase('./isabel64.png', true);
toc

tic
figGenSFCOctree3D('./data/SPH/db_4k_00025.csv');
toc

tic
figGenSFC2DCase('./OASIS/slice00198.csv');
toc
