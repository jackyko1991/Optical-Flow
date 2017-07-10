data_folder = dir('I:\optical flow\optical flow\data\30_fps');

for i = 10:10%(length(data_folder)-1)
    patient = data_folder(i).name
    path = 'I:\optical flow\optical flow\data\30_fps\',patient,'\';
    savepath = ['I:\optical flow\optical flow\data\30_fps\',patient,'\result_1\'];
    load([savepath 'tmp.mat']);
    
end

ROI_num = 3;    %determine number of velocity profile ROIs
ROI_disp_slice = 23;    %number of display slice
for i = 1:ROI_num;
    imshow(diff_0(:,:,ROI_disp_slice),[]);
    h = imrect(gca,[20 20 40 40]);
    ROI_position(i,:) = wait(h);
    clear h a
    close gcf
end

mask = zeros(size(diff_0));
for i = 65:66%size(image,3)
    %% image1, image2, alpha, iterations, initial Vx, initial Vy, display of flow (logic), dipslay of image
    disp(['Frame ' num2str(i) ' out of ' num2str(size(image,3)) ' is in progress'])
    
    %display (i-1)th image
    f1 = figure(1);
    imshow(image(:,:,i-1),[]);
    colormap(gray);
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f1,[savepath 'image_' num2str(i-1) '.jpg']);
    close all
    
    %display i-th image
    f2 = figure(2);
    imshow(image(:,:,i),[]);
    colormap(gray);
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f2,[savepath 'image_' num2str(i) '.jpg']);
    close all
    
    %display the subtraction images in subplot view
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
    
    % diplay overlayed magnitude image
    f4 = figure(4);
    mag = sqrt(Vx.^2+Vy.^2);
    [X, Y] = meshgrid(1:step:s(2), s(1):-step:1);
    Vx(mag > maxVel) = 0;
    Vy(mag > maxVel) = 0;
    %reduce display velocity display martix size
    u = interp2(Vx, X, Y);
    v = interp2(Vy, X, Y);
    mag_q = sqrt(u.^2+v.^2);
    image_overlay = RGB_part_overlay(image(:,:,i),mag,mask(:,:,i),maxVel);
    imshow(image_overlay);
    set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
    set(gca,'position',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(f4,[savepath 'mag_' num2str(i) '.jpg']);
    
    % display 
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
    
    % calculate velocity profile within ROI
    for j = 1:ROI_num
        ROI_crop = imcrop(image(:,:,i),ROI_position(j,:));
        mag_crop = imcrop(mag,ROI_position(j,:));
        mask_crop = logical(imcrop(mask(:,:,i),ROI_position(j,:)));
        avg(i,j) = mean(mag_crop(mask_crop));
    end
    
    close all;
    clear f1 f2 f3 f4 f5 f6
end

%% display ROI
cmap = hsv(size(avg,2)); %color map for plotting color
figure;
imshow(diff_0(:,:,ROI_disp_slice),[]);
% Enlarge figure to full screen.
set(gcf,'name','ROI selection','numbertitle','off');
set(gca,'position',[0 0 1 1]);
set(gcf,'PaperPositionMode','auto');
hold on
for i = 1:ROI_num
    h=imrect(gca,ROI_position(i,:));
    text(ROI_position(i,1)+ROI_position(i,3)/2,ROI_position(i,2)+ROI_position(i,4)/2,num2str(i),...
        'color',cmap(i,:));
end
saveas(gcf,[savepath 'ROI.jpg']);
close all

%% plot ROI velocity result
t = 0:1/30:(size(avg,1)-1)*1/30;
figure
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'name','Velocity Profile','numbertitle','off');
xlabel('Time(s)');
ylabel('Mean contrast flow velocity (pixel/s)');
hold on
for i = 1:size(avg,2)
    plot(t,avg(:,i).*30,'-o','Color',cmap(i,:));
    legendInfo{i} = ['ROI ' num2str(i)];
end
legend(legendInfo,'Location','southeast');

for i = 1:size(avg,2)
    hold on
    avg_tmp = avg(:,i).*30;
    avg_tmp(isnan(avg_tmp)) = 0;
    plot([t(1) t(end)],[mean(avg_tmp(avg_tmp~=0)) mean(avg_tmp(avg_tmp~=0))],'-.','Color',cmap(i,:),'LineWidth',2);
    disp(['Average value of ROI ' num2str(i) ' = ' num2str(mean(avg_tmp(avg_tmp~=0)))]);
    clear avg_tmp
end
saveas(gcf,[savepath 'ROI_velocity.jpg']);

clear t cmap
save([savepath 'tmp.mat']);