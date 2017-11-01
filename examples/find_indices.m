function [ ind ] = find_indices( address )
%% This script produces mask described in  
% Adapted Random Sampling Patterns for Accelerated MRI, 2011

%% Inputs:
% Acquire training data
% x_original = imread(address) % load your image 

x_original = address;
rate       = 1/5; % 20% of measurement.

% Output:
% the mask variable 'mask_pdf'

%% 

Nx              = size(x_original,1);
Ny              = size(x_original,2);
L               = Nx*Ny;
X_original      = abs(ifftshift(fft2(x_original)));
iPDF            = X_original./norm( X_original, 'fro' );

seq             = datasample(1:L,floor(rate*L)+1, 'Weights', iPDF(:), 'Replace', false);
mask_pdf        = zeros(Nx,Ny);
mask_pdf(seq)   = 1;

% imagesc(mask_pdf)
%sum(mask_pdf(:))/Nx/Ny;

ind             = find(mask_pdf);

end

