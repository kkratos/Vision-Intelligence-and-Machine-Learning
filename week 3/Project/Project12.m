% load the image we will experiment with
I = imresize(double(rgb2gray(imread('lena.png'))),[256 256]);

% build the Laplacian pyramid of this image with 6 levels
depth = 6;
L = laplacianpyr(I,depth);

% compute the quantization of the Laplacian pyramid
bins = [16,32,64,128,128,256]; % number of bins for each pyramid level
LC = encoding(L,bins);

% compute the entropy for the given quantization of the pyramid
ent = pyramident(LC);

function ent = pyramident(LC)

    % Input:
    % LC: the quantized version of the images stored in the Laplacian pyramid
    % Output:
    % br: the bitrate for the image given the quantization
    
    % Please follow the instructions to fill in the missing commands.
    
    ent = 0;                % initialization of entropy
    [r, c] = size(LC{1});
    pixI = r*c;             % number of pixels in the original image
    
    for i = 1:numel(LC)
        
        % 1) Compute the number of pixels at this level of the pyramid
        [rows, cols] = size(LC{i}); 
        pixel = rows * cols;
        % 2) Compute the entropy at this level of the pyramid 
        % (MATLAB command: entropy)
        ntropy = entropy(LC{i});
        % 3) Each level contributes to the entropy of the pyramid by a
        % factor that is equal to the sample density at this level, times
        % the entropy at this level. The sample density is computed as
        % (number of pixels at this level)/(number of pixels of original image).
        % Add this to the current sum of the entropy 'ent'
        ent = ent+(pixel/pixI)*ntropy;
    end
    
end

function LC = encoding(L, bins)

    % Input:
    % L: the Laplacian pyramid of the input image
    % bins: The number of bins used for discretization of each pyramid level
    % Output:
    % LC: the quantized version of the image stored in the Laplacian pyramid
    
    % Please follow the instructions to fill in the missing commands.

    depth = numel(bins);
    LC = cell(1,depth);
    
    for i = 1:depth

        % 1) Compute the edges of the bins we will use for discretization
        % (MATLAB command: linspace)
        % For level i, the linspace command will give you a row vector 
        % with bins(i) linearly spaced points between [X1,X2].
        % Remember that the range [X1,X2] depends on the level of the 
        % pyramid. The difference images (levels 1 to depth-1) are in 
        % the range of [-128,128], while the blurred image is in the 
        % range of [0,256]
        if i == depth % blurred image in range [0, 256]
            edges = linspace(0, 256, bins(i));
        else % difference image in range [-128,128]
            edges = linspace(-128, 128, bins(i));
        end
        
        % 2) Compute the centers that correspond to the above edges
        % The 1st center -> (1st edge + 2nd edge)/2
        % The 2nd center -> (2nd edge + 3rd edge)/2 and so on
        centers = (edges(1:end-1)+edges(2:end))/2;
        
        % 3) Discretize the values of the image at this level of the
        % pyramid according to edges (MATLAB command: discretize)
        % Hint: use 'centers' as the third argument of the discretize
        % command to get the value of each pixel instead of the bin index.
        LC{i} = discretize(L{i}, edges, centers)
        
    end
    
end

function I = collapse(L)

    % Input:
    % L: the Laplacian pyramid of an image
    % Output:
    % I: Recovered image from the Laplacian pyramid

    % Please follow the instructions to fill in the missing commands.
    
    depth = numel(L);
    
    % 1) Recover the image that is encoded in the Laplacian pyramid
    for i = depth:-1:1
        if i == depth
            % Initialization of I with the smallest scale of the pyramid
            I = L{i}
        else
            % The updated image I is the sum of the current level of the
            % pyramid, plus the expanded version of the current image I.
            I = expand(I) + L{i};
        end
    end

end

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