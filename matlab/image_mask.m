function image_mask

clear all
clc

% addpath LauSzeHungSamon_Pre;
addpath ChanTakYan_Post; % IMG-0003-00035.jpg
addpath ChanTakYan_Pre;  % IMG-0002-00030.jpg

MaxImage = [];
MinImage = [];

for k = 1:35
    %
    %
    if k <= 9
        Image = imread(['IMG-0003-0000',num2str(k),'.jpg']);
    else
        Image = imread(['IMG-0003-000',num2str(k), '.jpg']);
    end
    %
    if size(Image,3) ~= 1
        Image = rgb2gray(Image);
    end
    %
    Image = double(Image);
    %
    if k == 1
        MaxImage = zeros(size(Image));
        MinImage = zeros(size(Image)) + 1000;
    end
    %
    MaxImage = max(MaxImage,Image);
    MinImage = min(MinImage,Image);
end

figure(1); imagesc(MaxImage); axis equal; axis tight; colorbar
figure(2); imagesc(MinImage); axis equal; axis tight; colorbar
figure(3); imagesc(MaxImage - MinImage); axis equal; axis tight; colorbar

mask = zeros(size(MinImage));
DivImage = MaxImage - MinImage;
mask(DivImage > 11) = 1;
figure(4); imagesc(mask); axis equal; axis tight; colorbar




