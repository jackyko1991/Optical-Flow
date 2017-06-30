function data_analysis_ChanTakYan

clear all
clc
% % % 
% Chan Tak Yan
load ChanTakYan_Pre_Vx.mat;  PreVx  = Frame_Vx;
load ChanTakYan_Pre_Vy.mat;  PreVy  = Frame_Vy;
load ChanTakYan_Post_Vx.mat; PostVx = Frame_Vx;
load ChanTakYan_Post_Vy.mat; PostVy = Frame_Vy;

PreVxy  = sqrt( PreVx .^2 + PreVy .^2 );
PostVxy = sqrt( PostVx.^2 + PostVy.^2 );

MaxPreVxy = max(PreVxy,[],3); MaxPostVxy = max(PostVxy,[],3);
MinPreVxy = min(PreVxy,[],3); MinPostVxy = min(PostVxy,[],3);

MaxPreVxy  = (7.5 - 1)*(MaxPreVxy  - MinPreVxy);
MaxPostVxy = (7.5 - 1)*(MaxPostVxy - MinPostVxy);

%
% Noise Filtering
%

psf = fspecial('gaussian', 11, 3);
MaxPreVxy  = conv2(MaxPreVxy ,psf,'same');
MaxPostVxy = conv2(MaxPostVxy,psf,'same');
%
%

MinBar = min([min(MaxPreVxy(:)) min(MaxPostVxy(:))])-0.05;
MaxBar = max([max(MaxPreVxy(:)) max(MaxPostVxy(:))]);

MaskImage = imread('ChanTakYan_Post_Velocity_Edit.tif');
MaskImage = (MaskImage(:,:,[1:3]));
MaskImage = double(rgb2gray(MaskImage));
MaskImage(MaskImage >  1e-3) = 1;
MaskImage(MaskImage <= 1e-3) = 0;

Background = imread('ChanTakYan_Post.jpg');

figure(1),
imagesc(MaxPostVxy ,[MinBar MaxBar]); axis equal; axis tight; handle = colorbar, hold on
%
figure(1)
h = imagesc(Background), hold on
set(gca,'visible','off','color','none')
set(h,'AlphaData',1 - MaskImage)
%
set(get(handle,'Title'),'string','pixels/s','FontSize',12)
%
print('-depsc','-tiff','-r300','Velocity_ChanTakYan_Post')
print('-dtiff','-r200','Velocity_ChanTakYan_Post')
%
%
%
%
MaskImage = imread('ChanTakYan_Pre_Velocity_Edit.tif');
if size(MaskImage,3) == 4
    MaskImage = squeeze(MaskImage(:,:,[1:3]));
    MaskImage = double(rgb2gray(MaskImage));
else
    if size(MaskImage,3) == 3
        MaskImage = double(rgb2gray(MaskImage));
    else
        MaskImage = double(MaskImage);
    end
    
end
MaskImage(MaskImage >  1e-3) = 1;
MaskImage(MaskImage <= 1e-3) = 0;

Background = imread('ChanTakYan_Pre.jpg');
%
%
figure(2),
imagesc(MaxPreVxy,[MinBar MaxBar]); axis equal; axis tight; handle = colorbar, hold on
%
figure(2),
h = imshow(Background), hold on
set(gca,'visible','off','color','none')
set(h,'AlphaData',1 - MaskImage)
%
set(get(handle,'Title'),'string','pixels/s','FontSize',12)
% figure(5)
% imshow(MaskImage,[])
print('-depsc','-tiff','-r300','Velocity_ChanTakYan_Pre')
print('-dtiff','-r200','Velocity_ChanTakYan_Pre')



