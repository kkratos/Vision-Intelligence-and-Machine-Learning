% we load the two images we will blend
A = im2double(imread('orange.png'));
B = im2double(imread('apple.png'));

% mask that defines the blending region
R = zeros(512,512); R(:,257:512)=1;

% depth of the pyramids
depth = 5;

% 1) we build the Laplacian pyramids of the two images
LA = laplacianpyr(A,depth);
LB = laplacianpyr(B,depth);

% 2) we build the Gaussian pyramid of the selected region
GR = gausspyr(R,depth); 

% 3) we combine the two pyramids using the nodes of GR as weights
[LS] = combine(LA, LB, GR);


function [LS] = combine(LA, LB, GR)
    
    % Input:
    % LA: the Laplacian pyramid of the first image
    % LB: the Laplacian pyramid of the second image
    % GR: Gaussian pyramid of the selected region
    % Output:
    % LS: Combined Laplacian pyramid
    
    % Please follow the instructions to fill in the missing commands.
    
    depth = numel(LA);
    LS = cell(1,depth);    
    
    % 1) Combine the Laplacian pyramids of the two images.
    % For every level d, and every pixel (i,j) the output for the 
    % combined Laplacian pyramid is of the form:
    % LS(d,i,j) = GR(d,i,j)*LA(d,i,j) + (1-GR(d,i,j))*LB(d,i,j)
    for i = 1:depth
        % Put your code here
        LS{i} = GR{i}.*LA{i} + (1-GR{i}).*LB{i};
    end
end

function L = laplacianpyr(I,depth)

    % Input:
    % I: the input image
    % depth: number of levels of the Laplacian pyramid
    % Output:
    % L: a cell containing all the levels of the Laplacian pyramid
    
    % Please follow the instructions to fill in the missing commands.
    
    L = cell(1,depth);
    
    % 1) Create a Gaussian pyramid
    % Use the function you already created.
    G = gausspyr(I, depth);

    % 2) Create a pyramid, where each level is the corresponding level of
    % the Gaussian pyramid minus the expanded version of the next level of
    % the Gaussian pyramid.
    % Remember that the last level of the Laplacian pyramid is the same as
    % the last level of the Gaussian pyramid.
    for i = 1:depth
        if i < depth
            % same level of Gaussian pyramid minus the expanded version of next level
            L{i} = G{i} - expand(G{i+1});
        else
            % same level of Gaussian pyramid
            L{i} = G{i};
        end
    end
    
end

function G = gausspyr(I,depth)

    % Input:
    % I: the input image
    % depth: number of levels of the Gaussian pyramid
    % Output:
    % G: a cell containing all the levels of the Gaussian pyramid
    
    % Please follow the instructions to fill in the missing commands.
    
    G = cell(1,depth);
    
    % 1) Create a pyramid, where the first level is the original image
    % and every subsequent level is the reduced version of the previous level
    for i = 1:depth
        if i == 1
            % original image
            G{i} = I;
        else
            % reduced version of the previous level
            I = imresize(I, 0.5);
            G{i} = reduce(G{i-1});
            
        end
    end

end

function g = expand(I)

    % Input:
    % I: the input image
    % Output:
    % g: the image after the expand operation

    % Please follow the instructions to fill in the missing commands.
    
    % 1) Create the expanded image. 
    % The new image should be twice the size of the original image.
    % So, for an n x n image you will create an empty 2n x 2n image
    % Fill every second row and column with the rows and columns of the original image
    % i.e., 1st row of I -> 1st row of expanded image
    %       2nd row of I -> 3rd row of expanded image
    %       3rd row of I -> 5th row of expanded image, and so on
    [height, width, third] = size(I)
    new = zeros(2*height, 2*width, third);
    new(1:2:end, 1:2:end, 1:end) = I;
    % 2) Create a Gaussian kernel of size 5x5 and 
    % standard deviation equal to 1 (MATLAB command fspecial)
    gauss = fspecial('gaussian', [5 5], 1); 
    % 3) Convolve the input image with the filter kernel (MATLAB command imfilter)
    % Tip: Use the default settings of imfilter
    % Remember to multiply the output of the filtering with a factor of 4
    g = imfilter(new, gauss) .* 4;
end

function g = reduce(I)

    % Input:
    % I: the input image
    % Output:
    % g: the image after Gaussian blurring and subsampling

    % Please follow the instructions to fill in the missing commands.
    
    % 1) Create a Gaussian kernel of size 5x5 and 
    % standard deviation equal to 1 (MATLAB command fspecial)
    gaussian = fspecial('gaussian', [5 5], 1);
    % 2) Convolve the input image with the filter kernel (MATLAB command imfilter)
    % Tip: Use the default settings of imfilter
    g = imfilter(I, gaussian);
    % 3) Subsample the image by a factor of 2
    % i.e., keep only 1st, 3rd, 5th, .. rows and columns
    g = g(1:2:end, 1:2:end, 1:end);
end