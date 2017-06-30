% the programme is developed base on Horn-Schunck optical flow method
% Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
% 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337

%% prepare MATLAB workspace
clc;
close all;
clear all;

%% parameter settting (see Coarse2FineTwoFrames.m for the definition of the parameters)
step = 10; %vector plot step size
alpha = 0.012;
ratio = 0.75;
minWidth = 10;
maxVel = 30;%% maximum velocity for colormapping
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;
mu = 17; %% noise threshold value
crop = 0; %% logical decision for cropping
start = 1;
last = 202;
disp_slice = 120;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

%% loading data
data_folder = dir('I:\optical flow\data\30_fps\');
for i = 3:3%(length(data_folder)-1)
    patient = data_folder(i).name
    nii_folder = ['I:\optical flow\data\30_fps\',patient,'\nii\']
    savepath = ['I:\optical flow\data\30_fps\',patient,'\result\'];
    mkdir(savepath);
    nii = load_untouch_nii([nii_folder,'1.nii']);
    
    if crop == 1
        disp_img = rot90(nii.img(:,:,disp_slice));
        imshow(disp_img,[]);
        h = imrect;
        position = wait(h);

        count = 1;
        for i = start:last
            image(:,:,count) = imcrop(rot90(nii.img(:,:,i)),position);
            count = count + 1;
        end
        close(gcf);
        clear count h position
    
    elseif crop == 0
        start = 1;
        last = size(nii.img,3);
        image = rot90(nii.img(:,:,start:last));
    end
    save([savepath 'tmp.mat']);
end


% s = size(image);
% 
% diff = zeros(size(image));
% diff_itv = diff;
% diff_0 = diff; %difference with first image
% for i = 1:s(3)
%     diff_0(:,:,i) = image(:,:,i) - image(:,:,1);
% end

