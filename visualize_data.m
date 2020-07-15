clc
clear all
close all


%% Plot the velocity and the scalar at a single time point (make sure figure is full size for surf)

u = load('step10000_quick_u_data.txt');
v = load('step10000_quick_v_data.txt');
phi = load('step10000_quick_phi_data.txt');

figure
quiver(u,v)
xlim([0, 500])
ylim([0,100])

figure
surf(phi)
view(2)




%% plot u and v a 6 equally spaced time points to see the progression

u1 = load('step10000_quick_u_data.txt');
v1 = load('step10000_quick_v_data.txt');

u2 = load('step20000_quick_u_data.txt');
v2 = load('step20000_quick_v_data.txt');

u3 = load('step30000_quick_u_data.txt');
v3 = load('step30000_quick_v_data.txt');

u4 = load('step40000_quick_u_data.txt');
v4 = load('step40000_quick_v_data.txt');

u5 = load('step50000_quick_u_data.txt');
v5 = load('step50000_quick_v_data.txt');

u6 = load('step60000_quick_u_data.txt');
v6 = load('step60000_quick_v_data.txt');

xLeft = 0;
xRight = 500;
yBottom = 0;
yTop = 100;


% time progression of velocity
figure

subplot(2,3,1);
quiver(u1,v1)
title('time 10')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);

subplot(2,3,2); 
quiver(u2,v2)
title('time 20')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);

subplot(2,3,3); 
quiver(u3,v3)
title('time 30')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);

subplot(2,3,4); 
quiver(u4,v4)
title('time 40')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);

subplot(2,3,5); 
quiver(u5,v5)
title('time 50')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);

subplot(2,3,6); 
quiver(u6,v6)
title('time 60')
xlim([xLeft, xRight]);
ylim([yBottom, yTop]);



%% video of progression of velocity 

% save images

xLeft = 0;
xRight = 500;
yBottom = 0;
yTop = 100;


nums = [1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000, 27500, 30000, 32500, 35000, 37500, 40000, 42500, 45000, 47500, 50000, 52500, 55000, 57500, 60000];             

%%%%%% be sure to type in figure into command window then resize it


for i = 1:size(nums,2)
    
    file_num = nums(i);
    u = load(strcat('step', num2str(file_num), '_quick_u_data.txt'));
    v = load(strcat('step', num2str(file_num), '_quick_u_data.txt'));
    subplot(2,1,1)
    quiver(u,v);
    pbaspect([5 1 1])
    xlim([xLeft, xRight]);
    ylim([yBottom, yTop]);
    text(10, 90, strcat('time = ', num2str(0.001 * file_num)), 'FontSize', 15)
    title('quiver plot of velocity', 'FontSize', 15)
    
    
    phi = load(strcat('step', num2str(file_num), '_quick_phi_data.txt'));
    subplot(2,1,2)
    fig = surf(phi, 'edgecolor', 'none');
    caxis([0 1.2])
    view(2)
    pbaspect([5 1 1])
    xlim([xLeft, xRight]);
    ylim([yBottom, yTop]);
    title('surf plot of scalar being transported', 'FontSize', 15)
    
    
    F(i) = getframe(gcf);
    
    
    %%%%%% change folder location
    saveas(fig, strcat('images_lowres_quick_Re400/data', num2str(i), '.png'))
end

hold off
close all

%% create quiver video


writerObj = VideoWriter('myVideo1.avi');
writerObj.FrameRate = 1;
    
open(writerObj);

for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end

close(writerObj);
