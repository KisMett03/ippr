function [bw, bw_inverted, bw_noiseRemoved, bw_disk, bw_horizontal, bw_vertical, skeleton] = preprocess(inputImg)
    % Set default target height if not provided
    if nargin < 2
        targetHeight = 1200; % Standard height for processing
    end

    % Convert to grayscale if the image is in RGB
    if size(inputImg, 3) == 3
        img_gray = rgb2gray(inputImg);
    else
        img_gray = inputImg;
    end

    % Calculate scaling factor to maintain aspect ratio
    [originalHeight, ~] = size(img_gray);
    scaleFactor = targetHeight / originalHeight;

    % Resize the image using bicubic interpolation (better for handwriting)
    img_gray = imresize(img_gray, scaleFactor, 'bicubic');

    % Enhance contrast using CLAHE for better feature detection
    img_enhanced = adapthisteq(img_gray, 'ClipLimit', 0.02);

    % Apply sharpening to enhance edges
    img_enhanced = imsharpen(img_enhanced, 'Radius', 1.5, 'Amount', 0.6);

    % Step 2: Binarize using adaptive thresholding
    bw = imbinarize(img_enhanced, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.4);

    % Step 3: Invert the binary image
    bw_inverted = ~bw;

    % Step 4: Remove small objects (noise) - use a proportional threshold based on image size
    noiseThreshold = round(50 * scaleFactor);
    bw_noiseRemoved = bwareaopen(bw_inverted, noiseThreshold);

    % Step 5: Fill small holes before skeletonization to ensure continuous strokes
    bw_filled = imfill(bw_noiseRemoved, 'holes');  % Fill any holes in the image
    
    % Step 6: Apply disk-based opening for curved strokes preservation
    se_disk = strel('disk', max(1, round(1 * scaleFactor)));
    bw_disk = imopen(bw_filled, se_disk);

    % Step 7: Apply horizontal line-based closing for baseline features
    se_horizontal = strel('line', max(5, round(5 * scaleFactor)), 0);
    bw_horizontal = imclose(bw_filled, se_horizontal);

    % Step 8: Apply vertical line-based closing for letter connections
    se_vertical = strel('line', max(3, round(3 * scaleFactor)), 90);
    bw_vertical = imclose(bw_filled, se_vertical);

    % IMPROVED SKELETON EXTRACTION
    % Create skeleton after filling holes
    se_smooth = strel('disk', max(1, round(1 * scaleFactor)));
    bw_for_skeleton = imclose(bw_filled, se_smooth);  % Apply smoothing before skeletonization

    % Fill small holes again if necessary
    bw_for_skeleton = imfill(bw_for_skeleton, 8, 'holes');

    % Create skeleton with branch pruning
    minBranchLength = max(4, round(4 * scaleFactor));
    skeleton = bwskel(bw_for_skeleton, 'MinBranchLength', minBranchLength);

    % Optionally, apply additional cleaning or morphological operations to refine skeleton
    skeleton = bwmorph(skeleton, 'clean');  % Clean small branches or irregularities

    % Return all processed images
end
