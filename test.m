clear cam

cam = webcam(1);

CAM_WIDTH = 640;
CAM_HEIGHT = 480;
%preview(cam)


while true

img = snapshot(cam);
img = rgb2gray(img);

%C = fast9(im, 25, 100);
%C = corner(img);
%plot(CAM_WIDTH - C(:,1), CAM_HEIGHT - C(:,2),'r.');

C = detectFASTFeatures(img);

imshow(img); hold on;
plot(C.selectStrongest(50));


xlim([1 CAM_WIDTH])
ylim([1 CAM_HEIGHT])

end


