clear
clc
close all

files = dir(fullfile( '*.txt'));
fileNames = {files.name};
opening_height = str2double(extractAfter((extractBefore(fileNames,'mm')),'sym_'));
files(isnan(opening_height)) = [];
opening_height(isnan(opening_height)) = [];

[~,idx] = sort(opening_height);
%fileNames = fileNames{idx};
files = files(idx);
fileNames = {files.name};

n = length(fileNames);
for i=1:n
    origin_data(i).data = readmatrix(fileNames{i});
    origin_data(i).data(:,1) = [];%序号这一栏是不要的
    origin_data(i).data = origin_data(i).data(find(abs(origin_data(i).data(:,1))<=1),:);%进出口管道不重要，也不要
    pause(0.00001);
end

%% Sampling
dx = 0.001;
dy = 0.001;

o_x = origin_data(n).data(:,1);
o_y = origin_data(n).data(:,2);

nx = round((max(o_x)-min(o_x))/dx);
ny = round((max(o_y)-min(o_y))/dy);

xq = linspace(min(o_x),max(o_x),nx);
yq = linspace(min(o_y),max(o_y),ny);

[x_grid,y_grid] = ndgrid(xq,yq);
%%
for i =1:n
    tic
    data_pod(i).data = griddata(origin_data(i).data(:,1),origin_data(i).data(:,2),...
        origin_data(i).data(:,5),x_grid,y_grid,'linear');
    data_pod(i).data(isnan(data_pod(i).data))=0;
    pause(0.00001);
    toc
end
%% Result of sampling
for i =1:n
    pc = data_pod(i).data(:,:)';
    pc(pc==0) = nan;
    pcolor(squeeze(pc));
    % clim([1.3e6,2e6]);
    shading interp
    colormap;
    colorbar;
    pause(0.05);
end

%% POD
nxy = nx*ny;
ALL = zeros(n,nxy);
for i=1:1:n
    ALL(i,:) = reshape(data_pod(i).data,1,nxy);
end
U0x = mean(ALL,1);%均值
U_m = ALL-U0x;%去除均值

[U,S,phiU] = svd(U_m,'econ');%svd decomposition

An = U*S;%得到一个n阶的数
Ds = diag(S).^2/n;%特征值

%% Reconstruction
clc
nmodes = 4;%The first 4 modes are choosed
U_new = zeros(size(phiU));
for i =1:nmodes
    tic
    V{i}=An(:,i).*phiU(:,i)';
    U_new = U_new + V{i}';
    pause(0.0000001);
    toc
end
%% Reshape to x-y-z
for i =1:n
    tic
    U_n(i).p=reshape(U_new(:,i),nx,ny);
    U_n(i).p = U_n(i).p+reshape(U0x,nx,ny);
    pc = U_n(i).p;
    pc(pc==0) = nan;
    pcolor(x_grid,y_grid,squeeze(pc));
    %clim([1.3e6,2e6]);
    %clim([800,850]);
    shading interp
    colormap;
    colorbar;
    pause(0.0001);
    toc
end
%% calculate the error
error = mean(abs(([U_n.p]-[data_pod.data])./mean([U_n.p],'all')),'all');
%% Boundary of the valve
% clc
% boundary.wall = readmatrix('boundary_wall.txt');
% boundary.hole= readmatrix('boundary_hole.txt');
% boundary.all = [boundary.wall;boundary.hole];
% boundary_xy = boundary.all(find(abs(boundary.all(:,4))<=1e-3),:);
% boundary_xy = boundary_xy(find(abs(boundary_xy(:,2))<=1),:);
% bvalve = [boundary_xy(:,2),boundary_xy(:,3)];
% 
% hole_x = [-0.22,0.22,0.22,0.17,0.17,-0.17,-0.17,-0.22,-0.22];
% hole_y = [0.01,0.01,0.13,0.13,0.12,0.12,0.13,0.13,0.01];
% 
% in = inpolygon(boundary_xy(:,2),boundary_xy(:,3),hole_x,hole_y);
% boundary_hole = boundary_xy(in,:);
% boundary_valve = boundary_xy(~in,:);
% 
% up = find(boundary_valve(:,3)>0.12 | ((boundary_valve(:,2)>0) & (boundary_valve(:,3)>-0.1)));
% valve_up = boundary_valve(up,:);
% valve_down = setdiff(boundary_valve,valve_up,"rows");
% valve_up = sortrows(valve_up,2);
% valve_down = sortrows(valve_down,2);
% 
% valve_up_sorted = valve_up(1,:);
% remain_point = valve_up(2:end,:);
% while ~isempty(remain_point)
%     [idx,~] = knnsearch(remain_point(:,2:3),valve_up_sorted(end,2:3));
%     valve_up_sorted = [valve_up_sorted;remain_point(idx,:)];
%     remain_point(idx,:) = [];
% end
% 
% valve_down_sorted = valve_down(1,:);
% remain_point_down = valve_down(2:end,:);
% while ~isempty(remain_point_down)
%     [idx,~] = knnsearch(remain_point_down(:,2:3),valve_down_sorted(end,2:3));
%     valve_down_sorted = [valve_down_sorted;remain_point_down(idx,:)];
%     remain_point_down(idx,:) = [];
% end
% valve_sorted = [valve_up_sorted(5:end-10,:);flipud(valve_down_sorted(5:end-4,:));valve_up_sorted(5,:)];
% plot(valve_sorted(:,2),valve_sorted(:,3));
% %%
% in_mode = inpolygon(x_grid,y_grid,valve_sorted(:,2),valve_sorted(:,3));
% save('inmode.mat','in_mode');

read inmode.mat
%% Figure of modes
for i =1:4
    subplot(4,1,i)
    %pc = reshape(phiU(:,i),nx,ny);
    pc = reshape(V{i}(8,:),nx,ny);
    pc(~in_mode) = nan;
    pc = pc./1000000;

    pcolor(x_grid,y_grid,pc);
    shading interp;
    cmap = jet(1000);

    colormap(cmap);

    cb = colorbar;
    %cb.Title.String = '             Temperature / ℃';
    cb.Title.String = '             Pressure / MPa';

    set(gca, 'FontName', 'Times New Roman', 'FontSize',8);
    xlim([-1,1]);
    ylim([-0.4,0.3]);
    % pause(1);
    title(['(',char(96+i),')'],'Position',[-1.13,0.25]);
end
set(gcf, 'Units', 'centimeters', 'Position', [30,10,7,10]);
set(gcf,'Color',[1 1 1]);
%% Figure of flow field
clc
pc1 = U_n(5).p;
pc2 = data_pod(5).data;
pc = pc1-pc2;
pc1(~in_mode) = nan;
%pc = pc./1000000;
pc1(900:1100,100:250) = 0;
pcolor(x_grid,y_grid,pc1);


colormap(othercolor('BuDRd_18'));
cb = colorbar;
shading interp;
cb.Title.String = '             Velocity / m/s';
% xlim([-0.5,0.5]);
% ylim([-0.4,0.3]);
%clim([1.2,2]);
set(gca, 'FontName', 'Times New Roman', 'FontSize',8);
set(gcf,'Color',[1 1 1]);
set(gcf, 'Units', 'centimeters', 'Position', [30,10,7,4]);
