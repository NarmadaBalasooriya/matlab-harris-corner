%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SC556 - Assignment 03
% Author: Narmada Balasooriya
% Registration Number: PGIS/SC/MSC/CSC/17/06
% Harris corner detection implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rows, cols] = harris_corner(image)
    
    % In this harris_corner detection function sigma, radius and the
    % threshold values are predefined
    
    % input = image
    
    %% image details
    image_d = double(image);
    n_rows = size(image,1);
    n_cols = size(image,2);
    
    %% sobel operator for derivative of Gaussian
    dx = [-1 0 1; -2 0 2; -1 0 1];
    dy = [1 2 1; 0 0 0; -1 -2 -1];
    
    % alternative filter -> better than sobel operator
    dx_n = [-1 0 1; -1 0 1; -1 0 1];
    dy_n = dx_n';
    
    %% convolutional operation for Gaussian
    
    I_x = conv2(image_d, dx_n, 'same');
    I_y = conv2(image_d, dy_n, 'same');  
    
    %% derivative of Gaussian
    
    filter_size = 4;
    sigma = 0.5;
    G1 = fspecial('gaussian', [filter_size, filter_size], sigma);
    
    I_xx = imfilter(I_x.^2, G1);
    I_yy = imfilter(I_y.^2, G1);
    I_xy = imfilter(I_x.*I_y, G1); 
    
    %{
    can you filter2 as well
    
    I_xx = filter2(G1, I_x.^2);
    I_yy = filter2(G1, I_y.^2);
    I_xy = filter2(G1, I_x.*I_y);
    %}
    
    %% Computation of Mathematical expresssions
    k = 0.1;
    
    det_M = (I_xx.*I_yy - I_xy.^2);
    % this is alternative harris corner
    det_M_alt = (I_xx.*I_yy - I_xy.^2) - k*(I_xx + I_yy).^2;
    
    trace_M = I_xx + I_yy;
    
    % Compute the response of the detector at each pixel using alternative
    % det_M_alt
    R = det_M_alt - k*(trace_M).^2;
    
    
    %% Performing non-maximal suppression and threshold
    
    radius = 3; % an approximated value to consider the radius of region
    threshold = 0.05;
    mask_size = 2*radius;
    mask = ordfilt2(R, mask_size^2, ones(mask_size));
    
    % remove the points on the border
    remove_b = zeros(size(R));
    remove_b(radius:end-radius, radius:end-radius) = 1;
    
    % finding the maxima
    maximaa = (R==mask)&(R>=threshold)&remove_b;
        [rows, cols] = find(maximaa);
    
    figure(2);
    imshow(image); hold on;
    plot(rows, cols, '*r');
    
    
end