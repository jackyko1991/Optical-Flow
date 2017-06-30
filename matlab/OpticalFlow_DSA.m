function OpticalFlow_DSA

clear all
clc

% % % % % % % % addpath LauSzeHungSamon_Pre;
% % % % % % % addpath ChanTakYan_Post; % IMG-0003-00035.jpg
% % % % % % % addpath ChanTakYan_Pre;  % IMG-0002-00030.jpg
% % % % % % 
% % % % % % % addpath ChongShuFai_Post; % IMG-0005-00028.jpg
% % % % % % % addpath ChongShuFai_Pre;  % IMG-0004-00029.jpg
% % % % % % 
% % % % % % % addpath CheungChiFai_Post; % IMG-0008-00023.jpg
% % % % % % % addpath CheungChiFai_Pre; % IMG-0007-00021.jpg
% % % % % % 
% % % % % % addpath LawChungNgok_Post; % IMG-0010-00023.jpg
% % % % % % addpath LawChungNgok_Pre;  % IMG-0009-00022.jpg
% % % % % % 
% % % % % % alpha      = 5;
% % % % % % iterations = 15;
% % % % % % 
% % % % % % images = [];
% % % % % % 
% % % % % % % IMG-0001-00001.jpg
% % % % % % 
% % % % % % kmin = 1;
% % % % % % kmax = 23; % 29
% % % % % % 
% % % % % % Frame_Vx = [];
% % % % % % Frame_Vy = [];
% % % % % % 
% % % % % % for k = kmin : kmax-1
% % % % % %     %
% % % % % %     ncount = 1;
% % % % % %     %
% % % % % %     for j = k : k+1
% % % % % %         %
% % % % % %         if j <= 9
% % % % % %             Image = imread(['LawChungNgok_Post\IMG-0010-0000',num2str(j),'.jpg']);
% % % % % %         else
% % % % % %             Image = imread(['LawChungNgok_Post\IMG-0010-000',num2str(j),'.jpg']);
% % % % % %         end
% % % % % %         %
% % % % % %         if size(Image,3) ~= 1
% % % % % %             Image = rgb2gray(Image);
% % % % % %         end
% % % % % %         %
% % % % % %         %
% % % % % %         images(:,:,ncount) = double(Image);
% % % % % %         ncount = ncount + 1;
% % % % % %     end
% % % % % %     %
% % % % % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % % % % %     %
% % % % % %     uv(:,:,1) = Vx;
% % % % % %     uv(:,:,2) = Vy;
% % % % % %     %
% % % % % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % % % % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % % % % %     %
% % % % % %     %
% % % % % %     Frame_Vx(:,:,k) = Vx;
% % % % % %     Frame_Vy(:,:,k) = Vy;
% % % % % %     %
% % % % % % end
% % % % % % 
% % % % % % save LawChungNgok_Post_Vx.mat Frame_Vx
% % % % % % save LawChungNgok_Post_Vy.mat Frame_Vy
% % % % % % 
% % % % % % 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % images = [];
% % % % % % 
% % % % % % % IMG-0001-00001.jpg
% % % % % % 
% % % % % % kmin = 1;
% % % % % % kmax = 22;
% % % % % % 
% % % % % % Frame_Vx = [];
% % % % % % Frame_Vy = [];
% % % % % % 
% % % % % % for k = kmin : kmax-1
% % % % % %     %
% % % % % %     ncount = 1;
% % % % % %     %
% % % % % %     for j = k : k+1
% % % % % %         %
% % % % % %         if j <= 9
% % % % % %             Image = imread(['LawChungNgok_Pre\IMG-0009-0000',num2str(j),'.jpg']);
% % % % % %         else
% % % % % %             Image = imread(['LawChungNgok_Pre\IMG-0009-000',num2str(j),'.jpg']);
% % % % % %         end
% % % % % %         %
% % % % % %         if size(Image,3) ~= 1
% % % % % %             Image = rgb2gray(Image);
% % % % % %         end
% % % % % %         %
% % % % % %         %
% % % % % %         images(:,:,ncount) = double(Image);
% % % % % %         ncount = ncount + 1;
% % % % % %     end
% % % % % %     %
% % % % % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % % % % %     %
% % % % % %     uv(:,:,1) = Vx;
% % % % % %     uv(:,:,2) = Vy;
% % % % % %     %
% % % % % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % % % % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % % % % %     %
% % % % % %     %
% % % % % %     Frame_Vx(:,:,k) = Vx;
% % % % % %     Frame_Vy(:,:,k) = Vy;
% % % % % %     %
% % % % % % end
% % % % % % 
% % % % % % save LawChungNgok_Pre_Vx.mat Frame_Vx
% % % % % % save LawChungNgok_Pre_Vy.mat Frame_Vy
% % % % % % 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % addpath WongWaiMing_Post; % IMG-0006-00015.jpg
% % % % % % addpath WongWaiMing_Pre;  % IMG-0005-00011.jpg
% % % % % % 
% % % % % % alpha      = 5;
% % % % % % iterations = 15;
% % % % % % 
% % % % % % images = [];
% % % % % % 
% % % % % % % IMG-0001-00001.jpg
% % % % % % 
% % % % % % kmin = 1;
% % % % % % kmax = 15; % 29
% % % % % % 
% % % % % % Frame_Vx = [];
% % % % % % Frame_Vy = [];
% % % % % % 
% % % % % % for k = kmin : kmax-1
% % % % % %     %
% % % % % %     ncount = 1;
% % % % % %     %
% % % % % %     for j = k : k+1
% % % % % %         %
% % % % % %         if j <= 9
% % % % % %             Image = imread(['WongWaiMing_Post\IMG-0006-0000',num2str(j),'.jpg']);
% % % % % %         else
% % % % % %             Image = imread(['WongWaiMing_Post\IMG-0006-000',num2str(j),'.jpg']);
% % % % % %         end
% % % % % %         %
% % % % % %         if size(Image,3) ~= 1
% % % % % %             Image = rgb2gray(Image);
% % % % % %         end
% % % % % %         %
% % % % % %         %
% % % % % %         images(:,:,ncount) = double(Image);
% % % % % %         ncount = ncount + 1;
% % % % % %     end
% % % % % %     %
% % % % % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % % % % %     %
% % % % % %     uv(:,:,1) = Vx;
% % % % % %     uv(:,:,2) = Vy;
% % % % % %     %
% % % % % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % % % % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % % % % %     %
% % % % % %     %
% % % % % %     Frame_Vx(:,:,k) = Vx;
% % % % % %     Frame_Vy(:,:,k) = Vy;
% % % % % %     %
% % % % % % end
% % % % % % 
% % % % % % save WongWaiMing_Post_Vx.mat Frame_Vx
% % % % % % save WongWaiMing_Post_Vy.mat Frame_Vy
% % % % % % 
% % % % % % 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % images = [];
% % % % % % 
% % % % % % % IMG-0001-00001.jpg
% % % % % % 
% % % % % % kmin = 1;
% % % % % % kmax = 11;
% % % % % % 
% % % % % % Frame_Vx = [];
% % % % % % Frame_Vy = [];
% % % % % % 
% % % % % % for k = kmin : kmax-1
% % % % % %     %
% % % % % %     ncount = 1;
% % % % % %     %
% % % % % %     for j = k : k+1
% % % % % %         %
% % % % % %         if j <= 9
% % % % % %             Image = imread(['WongWaiMing_Pre\IMG-0005-0000',num2str(j),'.jpg']);
% % % % % %         else
% % % % % %             Image = imread(['WongWaiMing_Pre\IMG-0005-000',num2str(j),'.jpg']);
% % % % % %         end
% % % % % %         %
% % % % % %         if size(Image,3) ~= 1
% % % % % %             Image = rgb2gray(Image);
% % % % % %         end
% % % % % %         %
% % % % % %         %
% % % % % %         images(:,:,ncount) = double(Image);
% % % % % %         ncount = ncount + 1;
% % % % % %     end
% % % % % %     %
% % % % % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % % % % %     %
% % % % % %     uv(:,:,1) = Vx;
% % % % % %     uv(:,:,2) = Vy;
% % % % % %     %
% % % % % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % % % % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % % % % %     %
% % % % % %     %
% % % % % %     Frame_Vx(:,:,k) = Vx;
% % % % % %     Frame_Vy(:,:,k) = Vy;
% % % % % %     %
% % % % % % end
% % % % % % 
% % % % % % save WongWaiMing_Pre_Vx.mat Frame_Vx
% % % % % % save WongWaiMing_Pre_Vy.mat Frame_Vy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % addpath WongWingFung_Post; % IMG-0004-00021.jpg
% % % addpath WongWingFung_Pre;  % IMG-0003-00025.jpg
% % % 
% % % alpha      = 5;
% % % iterations = 15;
% % % 
% % % images = [];
% % % 
% % % % IMG-0001-00001.jpg
% % % 
% % % kmin = 1;
% % % kmax = 21; % 29
% % % 
% % % Frame_Vx = [];
% % % Frame_Vy = [];
% % % 
% % % for k = kmin : kmax-1
% % %     %
% % %     ncount = 1;
% % %     %
% % %     for j = k : k+1
% % %         %
% % %         if j <= 9
% % %             Image = imread(['WongWingFung_Post\IMG-0004-0000',num2str(j),'.jpg']);
% % %         else
% % %             Image = imread(['WongWingFung_Post\IMG-0004-000',num2str(j),'.jpg']);
% % %         end
% % %         %
% % %         if size(Image,3) ~= 1
% % %             Image = rgb2gray(Image);
% % %         end
% % %         %
% % %         %
% % %         images(:,:,ncount) = double(Image);
% % %         ncount = ncount + 1;
% % %     end
% % %     %
% % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % %     %
% % %     uv(:,:,1) = Vx;
% % %     uv(:,:,2) = Vy;
% % %     %
% % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % %     %
% % %     %
% % %     Frame_Vx(:,:,k) = Vx;
% % %     Frame_Vy(:,:,k) = Vy;
% % %     %
% % % end
% % % 
% % % save WongWingFung_Post_Vx.mat Frame_Vx
% % % save WongWingFung_Post_Vy.mat Frame_Vy
% % % 
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % images = [];
% % % 
% % % % IMG-0001-00001.jpg
% % % 
% % % kmin = 1;
% % % kmax = 25;
% % % 
% % % Frame_Vx = [];
% % % Frame_Vy = [];
% % % 
% % % for k = kmin : kmax-1
% % %     %
% % %     ncount = 1;
% % %     %
% % %     for j = k : k+1
% % %         %
% % %         if j <= 9
% % %             Image = imread(['WongWingFung_Pre\IMG-0003-0000',num2str(j),'.jpg']);
% % %         else
% % %             Image = imread(['WongWingFung_Pre\IMG-0003-000',num2str(j),'.jpg']);
% % %         end
% % %         %
% % %         if size(Image,3) ~= 1
% % %             Image = rgb2gray(Image);
% % %         end
% % %         %
% % %         %
% % %         images(:,:,ncount) = double(Image);
% % %         ncount = ncount + 1;
% % %     end
% % %     %
% % %     [Vx,Vy] = OpticalFlow(images,alpha,iterations);
% % %     %
% % %     uv(:,:,1) = Vx;
% % %     uv(:,:,2) = Vy;
% % %     %
% % %     figure(2*k - 1), plotflow(uv); title('Vector plot');
% % %     figure(2*k + 0), imagesc(sqrt(Vx.^2 + Vy.^2)); axis equal; axis tight; colorbar
% % %     %
% % %     %
% % %     Frame_Vx(:,:,k) = Vx;
% % %     Frame_Vy(:,:,k) = Vy;
% % %     %
% % % end
% % % 
% % % save WongWingFung_Pre_Vx.mat Frame_Vx
% % % save WongWingFung_Pre_Vy.mat Frame_Vy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath Post; % IMG-0004-00025.jpg
addpath Pre;  % IMG-0003-00020.jpg

alpha      = 5;
iterations = 15;

images = [];

% IMG-0001-00001.jpg

kmin = 1;
kmax = 25; % 29

Frame_Vx = [];
Frame_Vy = [];

uv = [];

for k = kmin : kmax-1
    %
    ncount = 1;
    %
    for j = k : k+1
        %
        if j <= 9
            Image = imread(['Post\IMG-0003-0000',num2str(j),'.jpg']);
        else
            Image = imread(['Post\IMG-0003-000',num2str(j),'.jpg']);
        end
        %
        if size(Image,3) ~= 1
            Image = rgb2gray(Image);
        end
        %
        %
        images(:,:,ncount) = double(Image);
        ncount = ncount + 1;
    end
    %
    [Vx,Vy] = OpticalFlow(images,alpha,iterations);
    %
    uv(:,:,1) = Vx;
    uv(:,:,2) = Vy;
    %
    f1 = figure;
    subplot(1,2,2)
    h = imshow(sqrt(Vx.^2 + Vy.^2),[]); 
    axis equal; 
    axis tight; 
    colormap(jet);
    %colorbar;
    subplot(1,2,1)
    plotflow(uv); 
    title('Vector plot');
    axis off;
    hold on;
    saveas(f1,['Post\IMG-0003-0000',num2str(j),'_result.jpg']);
    close all
    %
    %
    Frame_Vx(:,:,k) = Vx;
    Frame_Vy(:,:,k) = Vy;
    %
end

save Post_Vx.mat Frame_Vx
save Post_Vy.mat Frame_Vy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

images = [];

kmin = 1;
kmax = 20; % 29

Frame_Vx = [];
Frame_Vy = [];

for k = kmin : kmax-1
    %
    ncount = 1;
    %
    for j = k : k+1
        %
        if j <= 9
            Image = imread(['Pre\IMG-0002-0000',num2str(j),'.jpg']);
        else
            Image = imread(['Pre\IMG-0002-000',num2str(j),'.jpg']);
        end
        %
        if size(Image,3) ~= 1
            Image = rgb2gray(Image);
        end
        %
        %
        images(:,:,ncount) = double(Image);
        ncount = ncount + 1;
    end
    %
    [Vx,Vy] = OpticalFlow(images,alpha,iterations);
    %
    uv(:,:,1) = Vx;
    uv(:,:,2) = Vy;
    %
    f1 = figure;
    subplot(1,2,2)
    h = imshow(sqrt(Vx.^2 + Vy.^2),[]); 
    axis equal; 
    axis tight; 
    colormap(jet);
    %colorbar;
    subplot(1,2,1)
    plotflow(uv); 
    title('Vector plot');
    axis off;
    hold on;
    saveas(f1,['Post\IMG-0002-0000',num2str(j),'_result.jpg']);
    close all
    %
    %
    Frame_Vx(:,:,k) = Vx;
    Frame_Vy(:,:,k) = Vy;
    %
end

save Pre_Vx.mat Frame_Vx
save Pre_Vy.mat Frame_Vy



