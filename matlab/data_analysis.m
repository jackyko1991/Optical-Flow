function data_analysis

clear all
clc

% % Chan Siong Kee
% load ChanSiongKee_Pre_Vx.mat; PreVx = Frame_Vx(:,:,[4:end]);
% load ChanSiongKee_Pre_Vy.mat; PreVy = Frame_Vy(:,:,[4:end]);
% load ChanSiongKee_Post_Vx.mat; PostVx = Frame_Vx(:,:,[4:end]);
% load ChanSiongKee_Post_Vy.mat; PostVy = Frame_Vy(:,:,[4:end]);

% Chan Tak Yan
load CheungMayYuen_Pre_Vx.mat; PreVx = Frame_Vx(:,:,[2:18]);
load CheungMayYuen_Pre_Vy.mat; PreVy = Frame_Vy(:,:,[2:18]);
load CheungMayYuen_Post_Vx.mat; PostVx = Frame_Vx(:,:,[2:18]);
load CheungMayYuen_Post_Vy.mat; PostVy = Frame_Vy(:,:,[2:18]);

% % % Cheung Chi Fai ?????????
% % load CheungChiFai_Pre_Vx.mat; PreVx = Frame_Vx(:,:,[1:18]);
% % load CheungChiFai_Pre_Vy.mat; PreVy = Frame_Vy(:,:,[1:18]);
% % load CheungChiFai_Post_Vx.mat; PostVx = Frame_Vx(:,:,[1:18]);
% % load CheungChiFai_Post_Vy.mat; PostVy = Frame_Vy(:,:,[1:18]);

% % % Chong Shu Fai
% % load ChongShuFai_Pre_Vx.mat; PreVx = Frame_Vx(:,:,[4:end]);
% % load ChongShuFai_Pre_Vy.mat; PreVy = Frame_Vy(:,:,[4:end]);
% % load ChongShuFai_Post_Vx.mat; PostVx = Frame_Vx(:,:,[4:end]);
% % load ChongShuFai_Post_Vy.mat; PostVy = Frame_Vy(:,:,[4:end]);

% % % Lam Chung Kin
% % load LamChungKin_Pre_Vx.mat; PreVx = Frame_Vx;
% % load LamChungKin_Pre_Vy.mat; PreVy = Frame_Vy;
% % load LamChungKin_Post_Vx.mat; PostVx = Frame_Vx;
% % load LamChungKin_Post_Vy.mat; PostVy = Frame_Vy;

% % % Lau Sze Hung Samon
% % load LauSzeHungSamon_Pre_Vx.mat; PreVx = Frame_Vx;
% % load LauSzeHungSamon_Pre_Vy.mat; PreVy = Frame_Vy;
% % load LauSzeHungSamon_Post_Vx.mat; PostVx = Frame_Vx;
% % load LauSzeHungSamon_Post_Vy.mat; PostVy = Frame_Vy;

% % % % Law Chung Ngok
% % % load LawChungNgok_Pre_Vx.mat; PreVx = Frame_Vx;
% % % load LawChungNgok_Pre_Vy.mat; PreVy = Frame_Vy;
% % % load LawChungNgok_Post_Vx.mat; PostVx = Frame_Vx;
% % % load LawChungNgok_Post_Vy.mat; PostVy = Frame_Vy;

% % Wong Chi Keung
% load WongChiKeung_Pre_Vx.mat;PreVx = Frame_Vx(:,:,[4:end]);
% load WongChiKeung_Pre_Vy.mat;PreVy = Frame_Vy(:,:,[4:end]);
% load WongChiKeung_Post_Vx.mat; PostVx = Frame_Vx(:,:,[4:end]);
% load WongChiKeung_Post_Vy.mat; PostVy = Frame_Vy(:,:,[4:end]);

% % % Wong Yuk Lin
% % load WongYukLin_Pre_Vx.mat;   PreVx = Frame_Vx;
% % load WongYukLin_Pre_Vy.mat;   PreVy = Frame_Vy;
% % load WongYukLin_Post_Vx.mat; PostVx = Frame_Vx;
% % load WongYukLin_Post_Vy.mat; PostVy = Frame_Vy;

% % % % Wong Wing Fung
% % % load WongWingFung_Pre_Vx.mat;   PreVx = Frame_Vx;
% % % load WongWingFung_Pre_Vy.mat;   PreVy = Frame_Vy;
% % % load WongWingFung_Post_Vx.mat; PostVx = Frame_Vx;
% % % load WongWingFung_Post_Vy.mat; PostVy = Frame_Vy;

% % % Yiu Ki Kwong
% % load YiuKiKwong_Pre_Vx.mat; PreVx = Frame_Vx;
% % load YiuKiKwong_Pre_Vy.mat; PreVy = Frame_Vy;
% % load YiuKiKwong_Post_Vx.mat; PostVx = Frame_Vx;
% % load YiuKiKwong_Post_Vy.mat; PostVy = Frame_Vy;

PreVxy  = sqrt( PreVx .^2 + PreVy .^2 );
PostVxy = sqrt( PostVx.^2 + PostVy.^2 );

MaxPreVxy = max(PreVxy,[],3); MaxPostVxy = max(PostVxy,[],3);
MinPreVxy = min(PreVxy,[],3); MinPostVxy = min(PostVxy,[],3);

MaxPreVxy  = mean(PreVxy,3);%MaxPreVxy  - MinPreVxy;
MaxPostVxy = mean(PostVxy,3);%MaxPostVxy - MinPostVxy;

psf = fspecial('gaussian', 5, 1);
MaxPreVxy  = conv2(MaxPreVxy ,psf,'same');
MaxPostVxy = conv2(MaxPostVxy,psf,'same');

MinBar = min([min(MaxPreVxy(:)) min(MaxPostVxy(:))]);
MaxBar = max([max(MaxPreVxy(:)) max(MaxPostVxy(:))]);

figure(1),imagesc(MaxPreVxy ,[MinBar MaxBar]); axis tight; colorbar, title('Pre-Operation')
figure(2),imagesc(MaxPostVxy,[MinBar MaxBar]); axis tight; colorbar, title('Post-Operation')

MaxPreVxy  = 255 * MaxPreVxy./max(MaxPreVxy(:));
MaxPostVxy = 255 * MaxPostVxy./max(MaxPostVxy(:));

% imwrite(uint8(MaxPreVxy),'ChanTakYan_Pre_Velocity.tif')
% imwrite(uint8(MaxPostVxy),'ChanTakYan_Post_Velocity.tif')



