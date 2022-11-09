function [binarymask,gridx,gridy,gridz] = stl2mask(STL_filename)
%converts an .stl file to binary mask

%import the stl tree
tree = import_STL(STL_filename);

%predefine the grid points to force the voxels to be isotropic
%and the same number of voxels in each dimesion
mintree = min(tree.solidVertices{1});
maxtree = max(tree.solidVertices{1});

%anything much above this tends to give a memory error
nvox = 512;

%choose the voxel size based on the extent of the tree
voxsize = min((maxtree-mintree)/nvox);

gridx = linspace(mintree(1), mintree(1) + nvox*voxsize, nvox);
gridy = linspace(mintree(2), mintree(2) + nvox*voxsize, nvox);
gridz = linspace(mintree(3), mintree(3) + nvox*voxsize, nvox);

%change the format of tree so it can be passed to VOXELISE
tree.faces=tree.solidFaces{1};
tree.vertices=tree.solidVertices{1};

%voxelise the .stl to mask
binarymask = VOXELISE(gridx,gridy,gridz,tree);
%[binarymask,gridCOx,gridCOy,gridCOz] = VOXELISE(gridx,gridy,gridz,STL);