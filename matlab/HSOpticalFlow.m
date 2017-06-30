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
crop = 1; %% logical decision for cropping
start = 40;
last = 125;
disp_slice = 100;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

%% loading data
path = 'D:\optical flow\data\30_fps\01_ChanLamNuen\';
folder = dir([path 'nii\']);
savepath = [path 'result_2_crop\'];
mkdir(savepath);

%{
for i = 1:length(folder)-2
    image(:,:,i) = rgb2gray(imread([path '\data\' num2str(i) '.jpg']));
end
%}
nii = load_untouch_nii([path 'nii\2.nii']);

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
    clear count h position
    
elseif crop == 0
    image = rot90(nii.img(:,:,start:last));
end

s = size(image);

diff = zeros(size(image));
diff_itv = diff;
diff_0 = diff; %difference with first image
for i = 1:s(3)
    diff_0(:,:,i) = image(:,:,i) - image(:,:,1);
end

mask = zeros(size(diff_0));
for i = 2:size(image,3)
    %% image1, image2, alpha, iterations, initial Vx, initial Vy, display of flow (logic), dipslay of image
    disp(['Frame ' num2str(i) ' out of ' num2str(size(image,3)) ' is in progress'])
    
    f1 = figure(1);
    imshow(image(:,:,i-1),[]);
    colormap(gray);
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f1,[savepath 'image_' num2str(i-1) '.jpg']);
    close all
    
    f2 = figure(2);
    imshow(image(:,:,i),[]);
    colormap(gray);
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f2,[savepath 'image_' num2str(i) '.jpg']);
    close all
    
    f3 = figure(3);
    subplot(2,2,1)
    imshow(image(:,:,i-1),[]);
    colormap(gray);
    subplot(2,2,2)
    imshow(image(:,:,i),[]);
    colormap(gray);
    subplot(2,2,3)
    tmp = imcomplement(diff_0(:,:,i));
    tmp(tmp <= mu) = 0;
    tmp = medfilt2(tmp,[5 5]);
    imshow(tmp,[-50 50]);
    colormap(gray);
    mask(:,:,i) = logical(tmp);
    clear tmp
    subplot(2,2,4)
    diff(:,:,i) = image(:,:,i) - image(:,:,i-1);
    tmp = imcomplement(diff(:,:,i));
    tmp(abs(tmp) <= mu) = 0;
    tmp = medfilt2(tmp,[5 5]);
    imshow(tmp,[-50 50]);
    colormap(gray);
    clear tmp
    saveas(f3,[savepath 'diff_' num2str(i) '.jpg']);

    tic
    [Vx,Vy,warpI2] = Coarse2FineTwoFrames(image(:,:,i-1),image(:,:,i),para);
    toc
    uv(:,:,1) = Vx;
    uv(:,:,2) = Vy;
    %save([savepath 'result' num2str(i) '.mat'],'image','diff','Vx','Vy')
    
    f4 = figure(4);
    mag = sqrt(Vx.^2+Vy.^2);
    [X, Y] = meshgrid(1:step:s(2), s(1):-step:1);
    Vx(mag > maxVel) = 0;
    Vy(mag > maxVel) = 0;
    u = interp2(Vx, X, Y);
    v = interp2(Vy, X, Y);
    mag_q = sqrt(u.^2+v.^2);
    image_overlay = RGB_part_overlay(image(:,:,i),mag,mask(:,:,i),maxVel);
    imshow(image_overlay);
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f4,[savepath 'mag_' num2str(i) '.jpg']);
    
    f5 = figure(5);
    quiverwcolorbar(X, Y, u, v,3, 'bounds',[0 maxVel]);
    xlim([0 s(2)]);
    ylim([0 s(1)]);
    axis off;
    saveas(f5,[savepath 'vec_' num2str(i) '.jpg']);
    colorbar('off')
    hold on
    h5 = imshow(image(:,:,i),[]);
    set(h5,'alphadata',0.7);
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f5,[savepath 'vec_overlay_' num2str(i) '.jpg']);
    
    close all;
end

seg = max(mask,[],3);
seg = bwareaopen(seg, 50);
seg = imfill(seg,'holes');
seg = imclose(seg,strel('disk',5));
S=skeleton(seg);
BoundaryDistance = getBoundaryDistance(seg,0);
f1 = figure(1);
imshow(seg); hold on;
cmax = 0;
for i=1:length(S)
    X = [];
    Y = [];
    R = [];
    L=S{i};
    for j = 1:length(L)
        C{i}(j,1) = BoundaryDistance(round(L(j,1)),round(L(j,2)));
    end
    cmax = max(max(C{i}),cmax);
    cmap = jet(ceil(cmax));
    cplot = C{i}(:,1);
    figure(f1)
    scatter(L(:,2),L(:,1),3, cplot,'d');
    set(gca,'CLim',[0 ceil(cmax)]);
    colormap(cmap);
    clear cplot L
    %plot(L(:,2),L(:,1),'-','Color','k');
    
    %3D surface gerneration
    f2 = figure(2);
    X = S{i}(:,2);
    Y = S{i}(:,1);
    R = C{i}(:,1);
    [X_tube Y_tube Z_tube] = tubeplot(X,Y,zeros(size(X,1),1),R,R,40,[0 0 1]);
    subdivs = 40 + 1; 
    V=R*ones(1,subdivs);
    h = surf(X_tube,Y_tube,Z_tube,V,...
        'edgecolor','none');
    axis off;
    surf2stl([savepath 'surface' num2str(i) '.stl'],X_tube,Y_tube,Z_tube);
    close(f2);
    clear X_tube Y_tube Z_tube
    %{
    if i == 1
        X = S{i}(:,2);
        Y = S{i}(:,1);
        R = C{i}(:,1);
    else
        if X(end) == S{i}(start,2) && Y(end) == S{i}(start,1)
             X = [X; S{i}(:,2)];
             Y = [Y; S{i}(:,1)];
             R = [R; C{i}(:,1)];
        elseif X(end) == S{i}(end,2) && Y(end) == S{i}(end,1)           
            X = [X; flipud(S{i}(:,2))];
            Y = [Y; flipud(S{i}(:,1))];
            R = [R; flipud(C{i}(:,1))];
        end
    end
    %}
end

%pause
%close all
%clear i folder