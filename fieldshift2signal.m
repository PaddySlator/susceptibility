function [signal,signalall] = fieldshift2signal(fieldshiftmap, TE, ~)

GAMMA = 267.5e06;

%get the size of the field shift map
                                            
nvoxx = size(fieldshiftmap,1);
nvoxy = size(fieldshiftmap,2);
nvoxz = size(fieldshiftmap,3);
nTE = length(TE);

phaseshift = zeros(nvoxx, nvoxy, nvoxz, nTE);
signalall = zeros(nvoxx, nvoxy, nvoxz, nTE);

for t=1:nTE
    %calculate the phase shift
    phaseshift(:,:,:,t) = -GAMMA * fieldshiftmap * TE(t);
    %calculate the signal over all voxels in the field shift map
    signalall(:,:,:,t) = exp(1i * phaseshift(t));
end
                        
%average the signal in the imaging voxels
signal = squeeze(sum(signalall,[1 2 3]));
%get the magnitude of the signal



end