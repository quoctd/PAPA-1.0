function [ x_original, ind, m, n ] = get_dataset( data_option, ...
                                     small_data_size, imaddress, n_resize)

if nargin < 4, n_resize = []; end

if data_option == 1 % real image
    load('data/x_original_realdata.mat')
    load('data/indices_BrainPhantom.mat')
    
    % we have to resize the image because I have subsampling indices for a
    % 512x512 image
    % we have to normalize the image
    x_original      = double(x_original);
    x_original      = x_original/(norm(x_original(:)));
    
    % Prepare the subsampling indices
    ind             = find_indices(x_original);
    m               = size(x_original, 1);
    n               = size(x_original, 2);

elseif data_option == 2 % synthetic image
    load('data/BrainPhantom.mat')
    load('data/indices_BrainPhantom.mat')
    % Prepare the subsampling indices
    mask            = zeros(512);
    mask(ind)       = 1;
    mask            = fftshift(mask);
    ind             = find(mask);
    m               = size(x_original, 1);
    n               = size(x_original, 2);
    x_original      = double(x_original);
    
elseif data_option == 3 % small example
    m               = small_data_size(1); % the desired size for the phantom
    n               = small_data_size(2); % the desired size for the phantom
    
    x_original      = phantom(n);
    x_original      = double(x_original);
    x_original      = x_original/(norm(x_original(:)));
    x_original      = imresize(x_original, [m n]);
    ind = find_indices(x_original);

elseif data_option == 4 % use image provided by user
    im          = imread(imaddress);
    if length(size(im)) > 2, im_bw = rgb2gray(im); else im_bw = im; end
    if ~isempty(n_resize)
        im_bw   = imresize(im_bw, n_resize);
    end
    x_original  = double(im_bw);
    x_original  = x_original/(norm(x_original(:)));
    %x_original  = imresize(x_original, [200, 300]);
    [m, n]      = size(x_original);
    if m > 800 && n > 800
        x_original  = imresize(x_original, [650, 650]);
    end
    m           = size(x_original, 1);
    n           = size(x_original, 2);
    ind         = find_indices(x_original);
    fprintf('The size of image: [%d, %d].\n', m, n);
end

end

