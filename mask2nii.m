function mask2nii(ImBinary,x,y,z,filename)
% convert a mask (i.e. binary image) of a tree into a single nifti
% INPUTS: 
% mask - binary image with 1=inside tree, 0=outside
% x,y,z - grid points
%
%

%save the mask as nifti first to initialise the metadata
niftiwrite(ImBinary,filename,Compressed=true);

%get the metadata
maskinfo = niftiinfo(filename);

%get the voxel size from the grids
voxsize = [(abs(x(1)-x(end))/length(x)), ...
    (abs(y(1) - y(end))/length(y)), ...
    (abs(z(1) - z(end))/length(z))];

%update the voxelwise in the metadata
maskinfo.PixelDimensions = voxsize;

%resave the mask with the required metadata
niftiwrite(ImBinary,filename,maskinfo,Compressed=true);

end