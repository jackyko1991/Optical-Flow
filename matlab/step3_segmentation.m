path = 'D:\optical flow\data\30_fps\01_ChanLamNuen\';
savepath = [path 'result_2_crop\'];
load([savepath 'tmp.mat']);

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