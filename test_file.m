%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Narmada Balasooriya
% Harris corner detection implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the images -> Read the checkerboard image
image_t = checkerboard;

%% Matlab Harris corner features using computer vision toolbox
matlab_corners = detectHarrisFeatures(image_t);

%% my harris corner detection implementation
[rows, cols] = harris_corner(image_t);

%% plotting
figure(1);
imshow(image_t); hold on;
plot(matlab_corners.selectStrongest(50));
hold on;
plot(rows, cols, '*r');
%hold off;
