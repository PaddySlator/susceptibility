function [fieldshiftmap, binarymask] = mesh2fieldshift_FPM(filename, options)

% Estimate susceptibility-induced field shift map of mesh using 
% the finite perturber method of Pathak et al. 
%
% Note that the code forces the same number of voxels in each dimension - 
% which is necessary so that the FFT shift works properly
%
% Pathak et al., A novel technique for modeling susceptibility-based 
% contrast mechanisms for arbitrary microvascular geometries: 
% the finite perturber method, NeuroImage 2008, 
% 10.1016/j.neuroimage.2008.01.022 

% INPUTS:
% filename - either path to an STL file or path to nifti file
%
%
%
%
% OUTPUTs: 
% fieldshiftmap - the susceptibility-induced field shift map

disp('Finite Perturber Method for calculating susceptibility-induced field shift map')

%check if the input is a nifti (binary mask) or stl (mesh)
if endsWith(filename,'.nii.gz') || endsWith(filename,'.nii')
    filetype = 'nii';
elseif endsWith(filename,'.stl')
    filetype = 'stl';
else
    error('Incorrect file type! Input should be .stl, .nii, or .nii.gz')
end

%load and preprocess the input filename
if strcmp(filetype,'nii')
    %load the binary mask and metadata
    binarymask = niftiread(filename);
    binarymaskinfo = niftiinfo(filename);
    disp(['loaded the binary mask from nifti file: ' filename])

   
    %pad array to give the same number of voxels in each dimension
    nvox = max(size(binarymask));
    %add zeros at the end to pad
    binarymask = padarray(binarymask,nvox - size(binarymask),'post');
    disp('made the binary mask square by padding with zeros')

    %calculate the grid points - shift so they start at zero
    gridx = linspace(0, nvox*binarymaskinfo.PixelDimensions(1), nvox);
    gridy = linspace(0, nvox*binarymaskinfo.PixelDimensions(2), nvox);
    gridz = linspace(0, nvox*binarymaskinfo.PixelDimensions(3), nvox);    
elseif strcmp(filetype,'stl')
    [binarymask,gridx,gridy,gridz] = stl2mask(filename);
    disp(['loaded the tree from stl file: ' filename])
    disp(['converted tree to binary mask: ' filename])
end

%step 2 - calculate the 3D magnetic field from a SINGLE finite perturber
%(calculate at each grid point? seems like they assume that the single finite perturber is
%in the middle of the grid points?)

disp(['calculating susceptibility map for ' filename])

%strength of main magnetic field
B0 = 3;
%direction of main magnetic field
B0vec = [1 0 0];

%define the hematocrit fraction
HCT = 0.5;
%magnetic susceptibilty of deoxyhaemaglobin is 2.26 ppm
X_dHb = 2.26 * 10^-6;
%Fractional oxygenation
Y = 0.95;
%difference in magnetic susceptibility between vessel and tissue
deltaChi = HCT * (1 - Y) * X_dHb;

%what they use in the paper
%deltaChi = 3 * 10^-8;

disp('simulation values:')
disp(['B0 = ' num2str(B0) ' T (strength of main magnetic field)'])
disp(['B0vec = ' num2str(B0vec) ' (direction of main magnetic field)'])
disp(['hematocrit fraction = ' num2str(HCT)])
disp(['magnetic susceptibilty of deoxyhaemaglobin fixed at ' num2str(X_dHb) ' ppm'])
disp(['fractional oxygenation = ' num2str(Y)])
disp(['magnetic susceptibility difference between vessel and tissue = ' num2str(deltaChi) 'equation: HCT * (1 - Y) * X_dHb'])





%define the position of the single finite perturber (a cube the same size
%as the grid "voxels")
nvox = max(size(binarymask));
%get the voxel size from the grids
voxsize = [(abs(gridx(1)-gridx(end))/length(gridx)), ...
    (abs(gridy(1) - gridy(end))/length(gridy)), ...
    (abs(gridz(1) - gridz(end))/length(gridz))];

disp(['voxel size = ' num2str(voxsize)])

perturber_index = round(nvox/2);
pertuber_pos = [gridx(perturber_index) gridy(perturber_index) gridz(perturber_index)]; 
%"radius" of the perturber
a = mean(voxsize)/2;

disp(['"radius" of single finite perturber = ' num2str(a)])

%

nvoxx = length(gridx);
nvoxy = length(gridy);
nvoxz = length(gridz);

disp(['grid size = (' num2str(nvoxx) ',' num2str(nvoxy) ',' num2str(nvoxz) ')' ])



theta = zeros([nvoxx nvoxy nvoxz]);

deltaB = zeros([nvoxx nvoxy nvoxz]);
for x=1:nvoxx
    for y=1:nvoxy
        for z=1:nvoxz
            %get position of this grid point
            p = [gridx(x) gridy(y) gridz(z)];               
            %vector between this grid point and the single finite
            %perturber  
            rvec = p-pertuber_pos;
            %distance
            r = norm(rvec);
            %angle between rvec and main magnetic field
            theta(x,y,z) = atan2(norm(cross(rvec,B0vec)),dot(rvec,B0vec));
            
            deltaB(x,y,z) = (6/pi) * (deltaChi / 3) * (a^3 / r^3) * (3 * cos(theta(x,y,z))^2 - 1) * B0;
        end
    end
end

%susceptibility change is zero in the position of the single finite
%perturber
deltaB(perturber_index,perturber_index,perturber_index) = 0;
%normalised susceptibility
deltaBnorm = deltaB/(deltaChi*B0);
disp('done set up')



%step 3 - 3D FFT of the vascular structure (binarised image?)
fftbinmask = fftn(binarymask);
disp('calculated 3D FTT of the binary mask (i.e. vascular structure)')


%step 4 - 3D FFT of the single finite perturber field
%fftdeltaB = fftn(deltaB);
fftdeltaBnorm = fftn(deltaBnorm);
disp('calculated 3D FTT of the single finite perturber field')

%normalised version
%fftdeltaBnorm = fftn(deltaBnorm);


%step 5 - pointwise multiplication of the two 3D FFTs
%fftprod = fftbinmask .* fftdeltaB;
%normalised version
fftprodnorm = fftbinmask .* fftdeltaBnorm;

disp('calculated pointwise product of the two 3D FFTs')





%step 6 - inverse fft of this product
%need to fftshift this - not 100% sure why!
%fieldshiftmap = fftshift(ifftn(fftprod));
%normalised version
fieldshiftmap = fftshift(ifftn(fftprodnorm));

disp('calculated inverse FFT of this product')
disp('done and done!')













end



