%% Main function for feature extraction
function [group1Features, group2Features, group3Features] = feature_ext(bw, bw_inverted, bw_noiseRemoved, bw_disk, bw_horizontal, bw_vertical, skeleton, bw_combined)
    % Group 1: Cursive & Print Handwriting
    % Group 2: Block Letters & Slanted Handwriting
    % Group 3: Angular & Modern Calligraphy Handwriting
    % Input: Processed binary images (bw, bw_inverted, etc.)
    % Output: Structs of handwriting features for the 3 groups

    %% Group 1: Cursive and Print Handwriting Features

    % Cursive Features
    % Feature 1: Connected Components per Word (Letter Connectivity)
    group1Features.cursive.connectedComponents = extractConnectivity(bw_vertical);

    % Feature 2: Smooth Curvature (Few Sharp Corners)
    group1Features.cursive.smoothCurvature = extractSmoothCurvature(bw_disk);

    % Feature 3: Continuous Stroke (Few Pen Lifts / Endpoints)
    group1Features.cursive.continuousStroke = extractContinuousStroke(skeleton);

    % Print Features
    % Feature 1: Separate Letters (High Connected Component Count)
    group1Features.print.separateLetters = extractSeparateLetters(bw_noiseRemoved);

    % Feature 2: Upright Orientation (Minimal Slant)
    group1Features.print.uprightOrientation = extractMinimalSlant(bw_noiseRemoved);

    % Feature 3: Moderate Corners and Curves (Balanced Stroke Shapes)
    group1Features.print.balancedStrokeShapes = extractBalancedStrokeShapes(bw_combined);

    %% Group 2: Block Letters and Slanted Handwriting Features

    % Block Letters Features

    % Dominant Horizontal/Vertical Lines (strong indicator of block letters)
    group2Features.blockLetters.hvDominance = extractDominantHVLines(bw);

    % Extracts bounding boxes for individual letters and computes the mean
    group2Features.blockLetters.aspectConsistency = measureLetterAspectRatioConsistency(bw);

    % Measures consistency of letter aspect ratios based on the bounding boxes
    group2Features.blockLetters.heightConsistency = measureLetterHeightConsistency(bw);

    % Slanted Handwriting Features
    % Slant Angle Consistency (consistent tilt across writing)
    group2Features.slanted.avgSlantAngle = extractSlantAngle(bw);

    % Letter Tilt Uniformity (letters all tilted at similar angle)
    group2Features.slanted.letterTiltUniformity = extractLetterTiltUniformity(bw_inverted);

    % Vertical Stroke Count (fewer vertical strokes in slanted writing)
    group2Features.slanted.verticalStrokeCount = extractVerticalStrokeCount(bw_vertical);

    %% Group 3: Angular and Modern Calligraphy Handwriting Features

    % Angular Handwriting Features
    % Edge Orientation Variance (diverse stroke directions -> angular)
    group3Features.angular.edgeOrientationVariance = extractEdgeOrientationVariance(bw);

    % Low Circularity (angular letters have low roundness)
    group3Features.angular.circularity = extractCircularity(bw);

    % High Corner Density (many sharp corners in angular style)
    group3Features.angular.cornerDensity = extractCornerDensity(bw_combined);

    % Modern Calligraphy Features
    % Stroke Width Variation (thick-thin stroke contrasts)
    group3Features.calligraphy.strokeWidthVariation = extractStrokeWidthVariation(bw);

    % Flourishes (decorative strokes in calligraphy)
    group3Features.calligraphy.flourishes = extractFlourishes(bw_disk);

    % Smooth Curves (flowing, smooth curves in calligraphy)
    group3Features.calligraphy.smoothCurves = extractSmoothCurves(bw_disk);

    % Extended Letter accending and decending (calligraphy)
    group3Features.calligraphy.extendedAscDescMeasure = extractExtendedAscDesc(bw);
end

% Helper Functions for Segmentation

function wordSegments = segmentWords(bw)
    % segmentWords segments a binary image into word regions.
    % Uses horizontal dilation to connect letters into words.

    % Create a horizontal structuring element to connect nearby letters
    se = strel('line', 20, 0); % 20-pixel horizontal line (adjust as needed)
    bw_dilated = imdilate(bw, se);

    % Label connected components in the dilated image (these correspond to words)
    cc = bwconncomp(bw_dilated);
    stats = regionprops(cc, 'BoundingBox', 'Area');

    wordSegments = {};
    areaThreshold = 100; % ignore very small components (likely noise)

    for i = 1:length(stats)

        if stats(i).Area >= areaThreshold
            bbox = stats(i).BoundingBox;
            wordSegment = imcrop(bw, bbox);
            wordSegments{end + 1} = wordSegment;
        end

    end

end

function letterSegments = segmentHandwriting(bw)
    % segmentHandwriting segments a binary image into individual letter regions.
    % Uses connected component analysis to isolate letters (or connected letter groups).

    cc = bwconncomp(bw);
    stats = regionprops(cc, 'BoundingBox', 'Area');
    letterSegments = {};
    areaThreshold = 50; % minimum area to be considered a valid letter

    for i = 1:length(stats)

        if stats(i).Area >= areaThreshold
            bbox = stats(i).BoundingBox;
            letterSegment = imcrop(bw, bbox);
            letterSegments{end + 1} = letterSegment;
        end

    end

end

% Cursive Features Done
function connectivity = extractConnectivity(bw_vertical)
    % Compute connected components per word to gauge letter connectivity (cursive).
    wordSegments = segmentWords(bw_vertical);

    if isempty(wordSegments)
        % Fallback: if word segmentation fails, analyze whole image
        se = strel('disk', 3);
        bw_dilated = imdilate(bw_vertical, se);
        cc = bwconncomp(bw_dilated);
        rawConnectivity = cc.NumObjects;
    else
        % Calculate total connected components across words
        totalComponents = 0;
        numWords = length(wordSegments);

        for i = 1:numWords
            seg = wordSegments{i};
            % Dilate slightly to connect close strokes within the word
            se = strel('disk', 3);
            seg_dilated = imdilate(seg, se);
            cc = bwconncomp(seg_dilated);
            totalComponents = totalComponents + cc.NumObjects;
        end

        % Average connected components per word (lower for cursive writing)
        rawConnectivity = totalComponents / numWords;
    end

    % Normalize to [0,1] range
    % Lower values (closer to 1 component) indicate more cursive writing
    % Cap at 10 components per word for normalization
    maxComponents = 10;
    connectivity = max(0, min(1, (maxComponents - rawConnectivity) / (maxComponents - 1)));
end

function smoothCurvatureVal = extractSmoothCurvature(bw_disk)
    % Measures the smoothness of curves by comparing corner pixels to total boundary length.
    % A higher value indicates smoother, less angular curves.

    % Close small gaps to smooth out boundaries
    se = strel('disk', 2);
    bw_closed = imclose(bw_disk, se);

    % Optionally, apply slight Gaussian smoothing on the binary image
    bw_smoothed = imgaussfilt(double(bw_closed), 1) > 0.5;

    % Extract the boundary of the shapes
    boundaryImg = bwperim(bw_smoothed);

    % Use branchpoints on the boundary as a proxy for sharp corners
    cornerMask = bwmorph(boundaryImg, 'branchpoints');
    cornerCount = sum(cornerMask(:));

    % Compute total boundary length (number of boundary pixels)
    totalBoundaryPixels = sum(boundaryImg(:));

    if totalBoundaryPixels == 0
        smoothCurvatureVal = 0;
    else
        % Smoothness score = 1 - (corner density along the boundary)
        smoothCurvatureVal = 1 - (cornerCount / totalBoundaryPixels);
    end

end

function continuousStrokeVal = extractContinuousStroke(skeleton)
    % Assess continuity of strokes using the skeleton image.
    % Fewer endpoints relative to stroke length indicates more continuous strokes (cursive).

    % Remove small spurs from the skeleton to reduce noise
    skeleton_clean = bwmorph(skeleton, 'spur', 5); % remove spurs up to 5 pixels long

    % Optionally connect slightly broken strokes
    se = strel('disk', 1);
    skeleton_dilated = imdilate(skeleton_clean, se);

    % Count endpoints in the cleaned skeleton
    endpointsImg = bwmorph(skeleton_dilated, 'endpoints');
    endpointCount = sum(endpointsImg(:));
    skeletonLength = sum(skeleton_dilated(:));

    if skeletonLength == 0
        continuousStrokeVal = 0; % no strokes present
    else
        % Continuity score = 1 - (endpoints per skeleton length)
        continuousStrokeVal = 1 - (endpointCount / skeletonLength);
    end

end

% Print Features Done
function separateLetters = extractSeparateLetters(bw_noiseRemoved)
    % Measures how separated the letters are (print style has high separation).
    % We erode slightly to split connected letters, then count components.

    % Slight erosion to disconnect touching letters
    se_erode = strel('disk', 1);
    bw_eroded = imerode(bw_noiseRemoved, se_erode);

    % Reconstruct the image to original after erosion (keeps letters separated)
    bw_reconstructed = imreconstruct(bw_eroded, bw_noiseRemoved);

    % Segment into individual letters
    letterSegments = segmentHandwriting(bw_reconstructed);

    if isempty(letterSegments)
        separateLetters = 0;
        return;
    end

    % Determine valid components and capture each component's area
    validComponents = 0;
    areas = [];

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        stats = regionprops(seg, 'Area');

        if ~isempty(stats)
            segArea = max([stats.Area]);
            areas = [areas, segArea]; %#ok<AGROW>
        end

    end

    if isempty(areas)
        separateLetters = 0;
        return;
    end

    % Determine a minimum area threshold as a fraction of median letter size
    medianArea = median(areas);
    minArea = max(10, medianArea * 0.2);

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        stats = regionprops(seg, 'Area');

        if ~isempty(stats) && max([stats.Area]) > minArea
            validComponents = validComponents + 1;
        end

    end

    % Normalize by image size to account for scaling (larger images have more letters)
    totalPixels = numel(bw_noiseRemoved);
    rawValue = validComponents * 100 / sqrt(totalPixels);
    % Cap rawValue to a maximum of 10 to avoid skew and normalize to [0,1]
    separateLetters = min(rawValue, 10) / 10;
end

function uprightOrientation = extractMinimalSlant(bw_noiseRemoved)
    % Computes the average orientation of letters. Closer to 0 means more upright (print).
    letterSegments = segmentHandwriting(bw_noiseRemoved);
    orientations = [];

    if ~isempty(letterSegments)
        % Calculate orientation of each segment (letter or letter group)
        for i = 1:length(letterSegments)
            seg = letterSegments{i};
            stats = regionprops(seg, 'Orientation');

            if ~isempty(stats)
                % If multiple components in segment, take average orientation
                segOrientation = mean([stats.Orientation]);
                orientations = [orientations, segOrientation]; %#ok<AGROW>
            end

        end

        if ~isempty(orientations)
            rawOrientation = mean(orientations);
        else
            rawOrientation = 0;
        end

    else
        % Fallback: use entire image orientation if segmentation fails
        stats = regionprops(bw_noiseRemoved, 'Orientation');

        if numel(stats) > 1
            orientations = [stats.Orientation];
            rawOrientation = mean(orientations);
        elseif ~isempty(stats)
            rawOrientation = stats.Orientation;
        else
            rawOrientation = 0;
        end

    end

    % Normalize to [0,1] where 1 = perfectly upright, 0 = maximum slant (45 degrees or more)
    uprightOrientation = 1 - min(1, abs(rawOrientation) / 45);
end

function balancedStrokeShapes = extractBalancedStrokeShapes(bw_combined)
    % Evaluates if strokes have a balanced mix of curves and corners.
    % High value if neither sharp corners nor open-ended lines dominate.

    letterSegments = segmentHandwriting(bw_combined);

    if isempty(letterSegments)
        letterSegments = {bw_combined};
    end

    scores = zeros(length(letterSegments), 1);

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        % Compute the boundary of the segment
        boundaryImg = bwperim(seg);
        % Use boundary branchpoints as approximation for corners
        cornerMask = bwmorph(boundaryImg, 'branchpoints');
        cornerCount = sum(cornerMask(:));
        % Use boundary endpoints as a measure of stroke terminations (curviness)
        endpointMask = bwmorph(boundaryImg, 'endpoints');
        endpointCount = sum(endpointMask(:));
        totalBoundaryPixels = sum(boundaryImg(:));

        if totalBoundaryPixels == 0
            scores(i) = 0;
        else
            % A high score means moderate number of corners and endpoints
            scores(i) = 1 - (cornerCount / totalBoundaryPixels) - (endpointCount / totalBoundaryPixels);
        end

    end

    % Average across all letters for overall measure
    balancedStrokeShapes = mean(scores);
end

% Block Letters Features

function hvDominance = extractDominantHVLines(bw)
    % extractDominantHVLines measures how many edges/strokes are near 0 or 90 deg
    % indicative of block letters with strong horizontal/vertical lines.
    %
    % Input:
    %   bw - binary image of the handwriting (white on black)
    % Output:
    %   hvDominance - a ratio in [0..1], higher = more block-like HV lines

    % 1) Detect edges
    edgesImg = edge(bw, 'canny');

    % 2) Compute gradient orientation at edge pixels
    [~, Gdir] = imgradient(edgesImg);

    % Convert angles to [-180, 180]
    Gdir = mod(Gdir + 180, 360) - 180;

    % 3) Count total edge pixels
    totalEdgePixels = sum(edgesImg(:));

    if totalEdgePixels < 1
        hvDominance = 0;
        return;
    end

    angles = Gdir(edgesImg);

    % 4) Define a tolerance for "near horizontal" or "near vertical"
    tolerance = 15; % degrees

    % near horizontal => angle in [-tolerance, tolerance]
    isHorizontal = (angles >= -tolerance & angles <= tolerance);

    % near vertical => angle in [90-tolerance, 90+tolerance] or [-90-tolerance, -90+tolerance]
    isVertical = ((angles >= 90 - tolerance & angles <= 90 + tolerance) | ...
        (angles >= -90 - tolerance & angles <= -90 + tolerance));

    hvCount = sum(isHorizontal | isVertical);

    % 5) Ratio of HV edges to total edges
    hvDominance = hvCount / totalEdgePixels;
end

function aspectConsistency = measureLetterAspectRatioConsistency(bw)
    % Measures consistency of letter aspect ratios (width/height) for block letters.
    % - Uses segmentHandwriting to extract individual letter regions.
    % - For each letter, it picks the component with the largest area and computes
    %   its aspect ratio.
    % - It then computes the coefficient of variation (CV) of these ratios and
    %   converts it into a normalized consistency score in [0,1] (1 indicates perfect consistency).
    %
    % Input:
    %   bw - Binary image containing letters.
    % Output:
    %   aspectConsistency - Normalized consistency score in [0,1].

    letterSegments = segmentHandwriting(bw);

    if isempty(letterSegments)
        aspectConsistency = 1;
        return;
    end

    ratios = zeros(length(letterSegments), 1);

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        stats = regionprops(seg, 'BoundingBox', 'Area');

        if ~isempty(stats)
            % Choose the bounding box corresponding to the largest component
            areas = [stats.Area];
            [~, maxIdx] = max(areas);
            bbox = stats(maxIdx).BoundingBox; % [x, y, width, height]
            ratios(i) = bbox(3) / bbox(4); % Compute aspect ratio (width/height)
        end

    end

    meanRatio = mean(ratios);

    if meanRatio == 0
        aspectConsistency = 1;
    else
        % Lower coefficient of variation (CV) indicates higher consistency.
        cv = std(ratios) / meanRatio;
        aspectConsistency = max(0, min(1, 1 - cv));
    end

end

function heightConsistency = measureLetterHeightConsistency(bw)
    % Input:
    %   bw - Binary image containing letters.
    % Output:
    %   heightConsistency - Normalized consistency score in [0,1].
    %
    % This implementation applies morphological filtering to improve segmentation,
    % segments the image into individual letters, computes the height-to-width ratio
    % for each letter, removes outlier ratios, and then calculates the consistency as
    % 1 - (coefficient of variation).

    % Apply morphological filtering to improve segmentation
    se = strel('disk', 1);
    bw_filtered = imopen(bw, se);
    bw_filtered = imclose(bw_filtered, se);

    % Segment the image into individual letters using the filtered image
    letterSegments = segmentHandwriting(bw_filtered);

    if isempty(letterSegments)
        heightConsistency = 1; % No letters detected; assume perfect consistency.
        return;
    end

    % Compute height-to-width ratios for each valid letter
    ratios = [];
    
    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        stats = regionprops(seg, 'BoundingBox', 'Area');
        if ~isempty(stats)
            % Choose the component with the largest area
            areas = [stats.Area];
            [~, maxIdx] = max(areas);
            bbox = stats(maxIdx).BoundingBox; % [x, y, width, height]
            % Ignore letters with very small height
            if bbox(4) > 5 && bbox(3) > 0
                ratio = bbox(4) / bbox(3); % height-to-width ratio
                ratios = [ratios, ratio]; %#ok<AGROW>
            end
        end
    end

    if isempty(ratios)
        heightConsistency = 1;
        return;
    end

    % Remove outlier ratios using the IQR method
    medRatio = median(ratios);
    Q1 = quantile(ratios, 0.25);
    Q3 = quantile(ratios, 0.75);
    IQR = Q3 - Q1;
    lowerBound = Q1 - 1.5 * IQR;
    upperBound = Q3 + 1.5 * IQR;
    filteredRatios = ratios(ratios >= lowerBound & ratios <= upperBound);
    
    if isempty(filteredRatios)
        filteredRatios = ratios;
    end

    % Compute the coefficient of variation (CV) of the remaining ratios
    meanRatio = mean(filteredRatios);
    if meanRatio == 0
        heightConsistency = 1;
    else
        cv = std(filteredRatios) / meanRatio;
        % Convert CV to a consistency score (lower CV means higher consistency)
        heightConsistency = max(0, min(1, 1 - cv));
    end
end

% Slanted Handwriting Features Done
function avgSlantAngle = extractSlantAngle(bw)
    % Improved slant angle extraction using PCA on word segments.
    % Segment the image into words first.
    wordSegments = segmentWords(bw);

    if isempty(wordSegments)
        wordSegments = {bw}; % Fallback: use entire image
    end

    weightedAngles = [];
    segmentAreas = [];

    for i = 1:length(wordSegments)
        seg = wordSegments{i};
        % Get coordinates of foreground pixels (assuming text is white)
        [r, c] = find(seg);

        if numel(r) < 10
            continue;
        end

        % Form coordinate matrix [x y] where x = column, y = row
        coords = [c, r];
        % Compute PCA on the coordinates
        mu = mean(coords, 1);
        C = cov(double(coords));
        [V, ~] = eig(C);
        % Principal component is the eigenvector with maximum eigenvalue.
        [~, idx] = max(diag(C)); % This method is not reliable; use eigs or sort eigenvalues.
        % Instead, compute eigenvalues and sort them:
        [V, D] = eig(C);
        [eigVals, sortIdx] = sort(diag(D), 'descend');
        principalVec = V(:, sortIdx(1));

        % Compute angle (in degrees) between principal vector and horizontal axis
        angleRad = atan2(principalVec(2), principalVec(1));
        angleDeg = abs(rad2deg(angleRad));
        % For slanted handwriting, angles near 0 indicate near-horizontal; we expect higher values.
        % Weight this segment's angle by its number of foreground pixels.
        segArea = numel(r);
        weightedAngles(end + 1, 1) = angleDeg; %#ok<AGROW>
        segmentAreas(end + 1, 1) = segArea; %#ok<AGROW>
    end

    if isempty(weightedAngles)
        rawSlantAngle = 0;
    else
        % Compute area-weighted average of the slant angles.
        rawSlantAngle = sum(weightedAngles .* segmentAreas) / sum(segmentAreas);
    end

    % Normalize to [0,1] where 1 = maximum slant (90 degrees), 0 = no slant
    avgSlantAngle = min(1, rawSlantAngle / 90);
end

function tiltUniformity = extractLetterTiltUniformity(bw_inverted)
    % Measures how uniform the letter orientation is across the text.
    % High value if all letters have similar slant (consistent tilt).

    letterSegments = segmentHandwriting(bw_inverted);

    if isempty(letterSegments)
        % Fallback: consider whole image components if no segmentation
        stats = regionprops(bw_inverted, 'Orientation', 'Area');

        if isempty(stats)
            tiltUniformity = 1;
            return;
        end

        areas = [stats.Area];
        % Ignore tiny components below 5% of max area (noise or punctuation)
        minAreaThreshold = 0.05 * max(areas);
        validIdx = areas >= minAreaThreshold;
        if ~any(validIdx), validIdx = true(1, length(stats)); end
        orientations = [stats(validIdx).Orientation];
        % Standard deviation of orientation (lower = more uniform)
        stdOrient = std(orientations);
        tiltUniformity = max(0, 1 - (stdOrient / 45));
        return;
    end

    allOrientations = [];

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        stats = regionprops(seg, 'Orientation', 'Area');

        if ~isempty(stats)
            areas = [stats.Area];
            minAreaThreshold = 0.05 * max(areas);
            validIdx = areas >= minAreaThreshold;
            if ~any(validIdx), validIdx = true(1, length(stats)); end
            segOrientations = [stats(validIdx).Orientation];
            allOrientations = [allOrientations, segOrientations];
        end

    end

    if isempty(allOrientations)
        tiltUniformity = 1;
    else
        stdOrient = std(allOrientations);
        tiltUniformity = max(0, 1 - (stdOrient / 45));
    end

end

function verticalStrokeCount = extractVerticalStrokeCount(bw_vertical)
    % Counts tall vertical strokes (letters like l, t in upright text).
    % Slanted handwriting will have fewer purely vertical strokes.

    letterSegments = segmentHandwriting(bw_vertical);

    if isempty(letterSegments)
        % Fallback: count vertical strokes from whole image
        cc = bwconncomp(bw_vertical);
        stats = regionprops(cc, 'BoundingBox');
        count = 0;

        for i = 1:length(stats)
            bbox = stats(i).BoundingBox; % [x, y, width, height]

            if bbox(4) / (bbox(3) + eps) > 3
                % If height is much greater than width, treat as vertical stroke
                count = count + 1;
            end

        end

        rawCount = count;
    else
        totalCount = 0;

        for i = 1:length(letterSegments)
            seg = letterSegments{i};
            cc = bwconncomp(seg);
            stats = regionprops(cc, 'BoundingBox');
            count = 0;

            for j = 1:length(stats)
                bbox = stats(j).BoundingBox;

                if bbox(4) / (bbox(3) + eps) > 3
                    count = count + 1;
                end

            end

            totalCount = totalCount + count;
        end

        rawCount = totalCount;
    end

    % Normalize to [0,1] based on reasonable maximum count
    % Assuming 20 vertical strokes is a practical upper limit
    maxCount = 20;
    verticalStrokeCount = min(1, rawCount / maxCount);
end

% Angular Handwriting Features
function orientationStd = extractEdgeOrientationVariance(bw)
    % Measures variability of stroke edge orientations.
    % High standard deviation of edge directions indicates more angular, multi-directional strokes.

    letterSegments = segmentHandwriting(bw);

    if isempty(letterSegments)
        letterSegments = {bw};
    end

    stdValues = zeros(length(letterSegments), 1);

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        % Compute gradient magnitude and direction
        [Gmag, Gdir] = imgradient(seg);
        % Consider only strong edges for orientation analysis
        threshold = 0.1 * max(Gmag(:));
        idx = Gmag > threshold;
        orientations = Gdir(idx);
        orientations_rad = deg2rad(orientations);

        if isempty(orientations_rad)
            stdValues(i) = 0;
        else
            stdValues(i) = std(orientations_rad);
        end

    end

    % Average orientation variance across letters (radians)
    rawVariance = mean(stdValues);

    % Normalize to [0,1]
    % Theoretical maximum standard deviation of uniform angles is ~1.5 radians
    maxVariance = 1.5;
    orientationStd = min(1, rawVariance / maxVariance);
end

function circularity = extractCircularity(bw)
    % Calculates average circularity of components (1 = perfect circle, lower = more elongated/angular).
    % Angular style tends to have low circularity.

    % Morphologically smooth the image to reduce noise before measuring
    se = strel('disk', 2);
    bw_smooth = imopen(bw, se);
    bw_smooth = imclose(bw_smooth, se);

    letterSegments = segmentHandwriting(bw_smooth);

    if isempty(letterSegments)
        letterSegments = {bw_smooth};
    end

    allCirc = [];
    allAreas = [];

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        cc = bwconncomp(seg);

        if cc.NumObjects == 0
            continue;
        end

        stats = regionprops(cc, 'Area', 'Perimeter');
        areas = [stats.Area];
        perimeters = [stats.Perimeter];
        % Ignore very small components when calculating circularity
        minSize = 30;
        validIdx = areas > minSize & perimeters > 0;

        if ~any(validIdx)
            continue;
        end

        areas = areas(validIdx);
        perimeters = perimeters(validIdx);
        % Circularity = 4π * area / perimeter^2 (normalized to <= 1)
        individualCircularities = 4 * pi * areas ./ (perimeters .^ 2);
        individualCircularities = min(individualCircularities, 1);
        % Compute area-weighted circularity for this segment
        segCircularity = sum(individualCircularities .* areas) / sum(areas);
        allCirc = [allCirc, segCircularity];
        allAreas = [allAreas, sum(areas)];
    end

    if isempty(allCirc)
        circularity = 0;
    else
        % Overall circularity weighted by area of segments
        circularity = sum(allCirc .* allAreas) / sum(allAreas);
    end

end

function cornerDensity = extractCornerDensity(bw_combined)
    % Detects the density of corner points in the writing.
    % Higher corner density indicates a more angular style.

    letterSegments = segmentHandwriting(bw_combined);

    if isempty(letterSegments)
        letterSegments = {bw_combined};
    end

    allCornerDensity = zeros(length(letterSegments), 1);
    areas = zeros(length(letterSegments), 1);

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        % Smooth each segment to avoid counting noise as corners
        se = strel('disk', 2);
        bw_smooth = imopen(seg, se);
        bw_smooth = imclose(bw_smooth, se);
        % Use Harris/Shi-Tomasi corner detector on the segment
        corners = corner(bw_smooth, 'MinimumEigenvalue', 'QualityLevel', 0.1);
        totalCorners = size(corners, 1);
        % Area of the segment (for normalization)
        stats = regionprops(bw_smooth, 'Area');
        segArea = sum([stats.Area]);

        if segArea == 0
            allCornerDensity(i) = 0;
        else
            % Raw corner density = corners per sqrt(area) (to account for scale)
            cornerDensityRaw = totalCorners / sqrt(segArea);
            % Normalize: scale such that a dense corner count yields ~1
            allCornerDensity(i) = min(cornerDensityRaw / 0.05, 1);
        end

        areas(i) = segArea;
    end

    if sum(areas) == 0
        cornerDensity = 0;
    else
        % Weighted average by area
        cornerDensity = sum(allCornerDensity .* areas) / sum(areas);
    end

end

% Modern Calligraphy Features

function strokeWidthVariation = extractStrokeWidthVariation(bw)
    % Measures variation in stroke width across the text.
    % Higher value if there is a large difference between thick and thin strokes (as in calligraphy).

    segments = segmentHandwriting(bw);

    if isempty(segments)
        segments = {bw};
    end

    weightedVariationSum = 0;
    totalArea = 0;

    for i = 1:length(segments)
        seg = segments{i};
        % Compute area (foreground pixels) for weighting this segment
        areaSeg = sum(seg(:));
        totalArea = totalArea + areaSeg;
        % Use distance transform to estimate stroke thickness at each point
        distMap = bwdist(~seg);
        strokeWidths = distMap(distMap > 0);

        if isempty(strokeWidths)
            variation = 0;
        else
            meanWidth = mean(strokeWidths);
            stdWidth = std(strokeWidths);

            if meanWidth == 0
                variation = 0;
            else
                % Coefficient of variation of stroke width, scaled to [0,1]
                cv = stdWidth / meanWidth;
                variation = min(cv / 0.8, 1);
            end

        end

        weightedVariationSum = weightedVariationSum + variation * areaSeg;
    end

    if totalArea == 0
        strokeWidthVariation = 0;
    else
        % Area-weighted average variation
        strokeWidthVariation = weightedVariationSum / totalArea;
    end

end

function flourishes = extractFlourishes(bw_disk)
    % Detects decorative flourishes by looking for elongated or exaggerated strokes.
    % Uses skeleton branch points and component aspect ratios as indicators.

    segments = segmentHandwriting(bw_disk);

    if isempty(segments)
        segments = {bw_disk};
    end

    totalFlourishScore = 0;
    totalSegArea = 0;

    for i = 1:length(segments)
        seg = segments{i};
        % Apply morphological closing/opening to isolate prominent strokes
        se_close = strel('disk', 2);
        bw_closed = imclose(seg, se_close);
        se_open = strel('disk', 1);
        bw_processed = imopen(bw_closed, se_open);

        % Skeletonize to identify branches (complex intersections) and endpoints
        skel = bwskel(bw_processed);
        branchPoints = bwmorph(skel, 'branchpoints');
        endPoints = bwmorph(skel, 'endpoints');
        numBranchPoints = sum(branchPoints(:));
        numEndPoints = sum(endPoints(:));

        % Analyze connected components in the original segment
        cc = bwconncomp(seg);
        stats = regionprops(cc, 'BoundingBox', 'Area', 'Perimeter');

        % Compute median component height for reference
        heights = zeros(1, length(stats));

        for j = 1:length(stats)
            heights(j) = stats(j).BoundingBox(4);
        end

        heights = heights(heights > 5);

        if isempty(heights)
            segFlourish = 0;
        else
            medHeight = median(heights);
            flourishScore = 0;

            for j = 1:length(stats)
                area = stats(j).Area;
                perim = stats(j).Perimeter;
                height = stats(j).BoundingBox(4);

                if area < 20 % skip very small specks
                    continue;
                end

                % Elongation factor (eccentricity proxy): higher if component is long/thin
                elongation = (perim ^ 2) / (4 * pi * area);
                % Height relative to median: captures unusually tall strokes (e.g., extended ascenders)
                heightRatio = height / medHeight;
                % If a component is highly elongated or much taller than typical, count it as flourish element
                if elongation > 3 || heightRatio > 1.5
                    flourishScore = flourishScore + area * (elongation / 10 + heightRatio / 2);
                end

            end

            % Normalize flourish score by segment area
            totalAreaSeg = numel(seg);
            normalizedScore = min(flourishScore / (totalAreaSeg * 0.05), 1);
            % Incorporate branchpoint-to-endpoint ratio as part of flourish indicator
            branchToEndRatio = numBranchPoints / max(1, numEndPoints);
            segFlourish = normalizedScore * 0.7 + min(branchToEndRatio / 3, 1) * 0.3;
            segFlourish = min(segFlourish, 1);
        end

        % Weight this segment’s flourish score by its size
        segArea = sum(seg(:));
        totalFlourishScore = totalFlourishScore + segFlourish * segArea;
        totalSegArea = totalSegArea + segArea;
    end

    if totalSegArea == 0
        flourishes = 0;
    else
        flourishes = totalFlourishScore / totalSegArea;
    end

end

function smoothCurves = extractSmoothCurves(bw_disk)
    % Measures how smooth and flowing the curves are in the text.
    % High value if boundaries of letters have low curvature variation (smooth strokes).

    segments = segmentHandwriting(bw_disk);

    if isempty(segments)
        segments = {bw_disk};
    end

    totalSmoothScore = 0;
    totalArea = 0;

    for i = 1:length(segments)
        seg = segments{i};
        % Find boundaries of the segment
        boundaries = bwboundaries(seg);
        segTotalCurvature = 0;
        segTotalPoints = 0;

        for k = 1:length(boundaries)
            boundary = boundaries{k};

            if size(boundary, 1) < 10
                continue; % skip very small boundaries
            end

            x = boundary(:, 2);
            y = boundary(:, 1);
            % Smooth the boundary coordinates to reduce pixel noise
            xSmooth = smooth(x, 5);
            ySmooth = smooth(y, 5);
            n = length(xSmooth);
            curvature = zeros(n, 1);
            % Compute curvature (angle change) along the boundary
            for j = 2:(n - 1)
                dx1 = xSmooth(j) - xSmooth(j - 1);
                dy1 = ySmooth(j) - ySmooth(j - 1);
                dx2 = xSmooth(j + 1) - xSmooth(j);
                dy2 = ySmooth(j + 1) - ySmooth(j);
                angle1 = atan2(dy1, dx1);
                angle2 = atan2(dy2, dx2);
                angleDiff = abs(angle2 - angle1);
                % Normalize angle difference to [0, π]
                angleDiff = mod(angleDiff + pi, 2 * pi) - pi;
                curvature(j) = abs(angleDiff);
            end

            % Handle endpoints of the boundary array
            curvature(1) = curvature(2);
            curvature(n) = curvature(n - 1);
            segTotalCurvature = segTotalCurvature + sum(curvature);
            segTotalPoints = segTotalPoints + n;
        end

        if segTotalPoints == 0
            segSmoothScore = 0;
        else
            % Average curvature per point; map lower curvature to higher smoothness score
            avgCurvature = segTotalCurvature / segTotalPoints;
            segSmoothScore = max(0, 1 - (avgCurvature / (pi / 2)));
        end

        % Weight by segment size
        areaSeg = sum(seg(:));
        totalSmoothScore = totalSmoothScore + segSmoothScore * areaSeg;
        totalArea = totalArea + areaSeg;
    end

    if totalArea == 0
        smoothCurves = 0;
    else
        smoothCurves = totalSmoothScore / totalArea;
    end

end

function extendedAscDescMeasure = extractExtendedAscDesc(bw)
    % extractExtendedAscDesc measures how much letters extend beyond
    % a core "x-height", typical of modern calligraphy ascenders & descenders.
    %
    % Input:
    %   bw - Binary image (white text on black) of the handwriting
    %
    % Output:
    %   extendedAscDescMeasure - [0..1] measure of extended ascenders/descenders

    letterSegments = segmentHandwriting(bw); % from your existing code

    if isempty(letterSegments)
        extendedAscDescMeasure = 0;
        return;
    end

    ratioSum = 0;
    totalLetters = 0;

    for i = 1:length(letterSegments)
        seg = letterSegments{i};
        % bounding box
        stats = regionprops(seg, 'BoundingBox', 'Area');

        if isempty(stats)
            continue;
        end

        % If multiple components in this segment, pick the largest area
        [~, maxIdx] = max([stats.Area]);
        bbox = stats(maxIdx).BoundingBox; % [x, y, width, height]
        letterHeight = bbox(4);

        % Heuristic: approximate x-height by ignoring top/bottom ~10%
        % so "core" height ~ 80% of bounding box
        % The ratio: extended portion / core portion
        % If letterHeight is small or near 0, skip.
        if letterHeight < 5
            continue;
        end

        coreHeight = 0.8 * letterHeight;
        extendedPortion = letterHeight - coreHeight; % top + bottom
        % ratio > 0 means some ascender/descender
        ratio = extendedPortion / coreHeight; % e.g. 0.5 means ascenders are 50 % of x-height
        ratioSum = ratioSum + ratio;
        totalLetters = totalLetters + 1;
    end

    if totalLetters == 0
        rawRatio = 0;
    else
        rawRatio = ratioSum / totalLetters;
    end

    % Normalize to [0,1]
    % A ratio of 1.5 is considered very extended
    maxRatio = 1.5;
    extendedAscDescMeasure = min(1, rawRatio / maxRatio);
end
