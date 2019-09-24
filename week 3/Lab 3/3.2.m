% loading the image
A = im2double(imread('orange.png'));
% depth of the pyramids
depth = 5;

% we build the Gaussian pyramid
GA = gausspyr(A,depth);

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

function g = reduce(I)

    % Add your code from the previous step
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