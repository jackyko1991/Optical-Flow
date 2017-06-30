data_folder = dir('I:\optical flow\data\30_fps\');

for i = 3:3%(length(data_folder)-1)
    patient = data_folder(i).name
    path = 'I:\optical flow\data\30_fps\01_ChanLamNuen\';
    savepath = ['I:\optical flow\data\30_fps\',patient,'\result\'];
    load([savepath 'tmp.mat']);
    
    for j = 65:65%size(image,3)
        tic
        [Vx,Vy,warpI2] = Coarse2FineTwoFrames(image(:,:,i-1),image(:,:,i),para);
        toc
        uv(:,:,1) = Vx;
        uv(:,:,2) = Vy;
        %save([savepath 'result' num2str(i) '.mat'],'image','diff','Vx','Vy')

        % diplay overlayed magnitude image
        f = figure;
        mag = sqrt(Vx.^2+Vy.^2);

        imshow(mag);
%         set(gcf,'position',[pos(1), pos(2), s(2), s(1)]);
        set(gca,'position',[0 0 1 1]);
        set(gcf,'PaperPositionMode','auto');
%         saveas(f4,[savepath 'mag_' num2str(i) '.jpg']);
        %close all
    end
end