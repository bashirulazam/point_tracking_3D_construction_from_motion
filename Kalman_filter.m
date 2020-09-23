clear all
close all
clc

fileprefix = 'proj4_img_seq/';
L = 9;
frame0 = double(imread(strcat(fileprefix,num2str(10),'.pgm')));
[x1,y1,vx1,vy1,x0,y0,vx0,vy0,bound_box] = get_init_cor(frame0,L);
imshow(uint8(frame0));
for i = 1:29
    im_seq(:,:,i) = double(imread(strcat(fileprefix,num2str(i+10),'.pgm')));       
end

Sigma(:,:,1) = [100 0 0 0;
                0 100 0 0;
                0  0 25 0;
                0  0  0 25;];
s(:,:,1) = [x1 y1 vx1 vy1]';            
for i = 1:28
        nextFrame = im_seq(:,:,i+1);
        for p = 1:length(x1)
            [s(:,p,i+1),Sigma(:,:,i+1)] = MyUpdate(s(:,p,i),Sigma(:,:,i),nextFrame,bound_box(:,:,p),L,frame0);
            x = round(s(1,p,i+1));
            y = round(s(2,p,i+1));
            im_seq(y,x,i+1) = 255;
        end
        figure 
        imshow(uint8(im_seq(:,:,i+1)));
        filename = strcat('Results/',num2str(i),'.jpg');
        imwrite(uint8(im_seq(:,:,i+1)),filename);
        
end

coords = zeros(2,size(s,2),size(s,3)+1);
coords(:,:,2:30) = round(s(1:2,:,:));
coords(:,:,1) = [x0'; y0'];
savefile = 'coords.mat';
save(savefile,'coords');
figure
for i = 1:29 
    Sigma_trace(i) = trace(Sigma(:,:,i))
end
plot(Sigma_trace)
xlabel('Frame Number');
ylabel('Trace of Covariance Matrices');
title('Plot of trace of Covariance Matrices');
axis([1 29 -10 250])
