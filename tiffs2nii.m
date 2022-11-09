function tiffs2nii(tiff_basename)
%convert a stack of tiff files into a single nifti
%tiff files should be numbered with leading zeros
% examples: filename0001.tiff to filename1000.tiff,
% filename001.tiff to filename100.tiff

tiff_files = dir([tiff_basename '*tiff']);

n_img = length(tiff_files);
%load the first tiff to get dimensions
thistiff = imread(tiff_files(1).name);

tiff_stack = zeros([size(thistiff) n_img]);

for i=2:n_img
    this_tiff = imread(tiff_files(i).name);
    tiff_stack(:,:,i) = this_tiff;
end

niftiwrite(tiff_stack, tiff_basename , Compressed=true);

end