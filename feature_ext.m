%% Main function for feature extraction
function [group1Features, group2Features, group3Features] = feature_ext(bw, bw_inverted, bw_noiseRemoved, bw_disk, bw_horizontal, bw_vertical, skeleton, bw_combined)

    % Group 1: Cursive & Print Handwriting
    % Group 2: Block Letters & Slanted Handwriting
    % Group 3: Angular & Modern Calligraphy Handwriting
    % Input: Processed binary images (bw, bw_inverted, etc.)
    % Output: Handwriting features for the 3 groups

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
    % Horizontal/Vertical Line Presence (Block letters are typically made up of straight lines)
    group2Features.blockLetters.linePresence = extractLinePresence(bw_vertical, bw_horizontal);

    % Right-Angle Corners (Block letters tend to have many sharp right angles)
    group2Features.blockLetters.angleCorners = extractRightAngleCorners(bw_combined);

    % Uniform Stroke Lengths (Block letters often have consistent stroke lengths)
    group2Features.blockLetters.strokeLengthConsistency = extractStrokeLengthConsistency(bw);

    % Slanted Handwriting Features
    % Slant Angle Consistency (Slanted handwriting tends to have consistent tilt)
    group2Features.slanted.avgSlantAngle = extractSlantAngle(bw);

    % Letter Tilt Uniformity (Slanted handwriting will have consistent letter tilt)
    group2Features.slanted.letterTiltUniformity = extractLetterTiltUniformity(bw_inverted);

    % Vertical Stroke Detection (Slanted handwriting has fewer vertical strokes)
    group2Features.slanted.verticalStrokeCount = extractVerticalStrokeCount(bw_vertical);

    %% Group 3: Angular and Modern Calligraphy Handwriting Features

    % Angular Handwriting Features
    % Corner Detection (Angular handwriting has sharp corners and acute angles)
    group3Features.angular.cornerCount = extractCornerCount(bw_combined);

    % Low Circularity (Angular letters have lower circularity due to angular shapes)
    group3Features.angular.circularity = extractCircularity(bw);

    % High Corner Density (Angular handwriting has many sharp corners in each letter)
    group3Features.angular.cornerDensity = extractCornerDensity(bw_combined);

    % Modern Calligraphy Features
    % Stroke Width Variation (Modern calligraphy has thick-thin stroke contrasts)
    group3Features.calligraphy.strokeWidthVariation = extractStrokeWidthVariation(bw);

    % Flourishes (Calligraphy often includes decorative flourishes)
    group3Features.calligraphy.flourishes = extractFlourishes(bw_disk);

    % Smooth Curves (Modern calligraphy has smooth, flowing curves)
    group3Features.calligraphy.smoothCurves = extractSmoothCurves(bw_disk);

end


function connectedComponents = extractConnectivity(bw_vertical)
    % Label connected components
    cc = bwconncomp(bw_vertical);
    
    % Get stats
    stats = regionprops(cc, 'Area');
    
    % Filter out small components
    minArea = 10;
    validComponents = [stats.Area] > minArea;
    
    % Get number of valid components
    connectedComponents = sum(validComponents);
    
    % Normalize by image area
    totalPixels = numel(bw_vertical);
    connectedComponents = connectedComponents / totalPixels;
end

function smoothCurvatureVal = extractSmoothCurvature(bw_disk)
    % Extract boundaries
    [B, L] = bwboundaries(bw_disk, 'noholes');
    
    % Initialize
    totalCurvature = 0;
    validBoundaryCount = 0;
    
    % Loop through each boundary
    for k = 1:length(B)
        boundary = B{k};
        
        % Skip small boundaries
        if size(boundary, 1) < 20
            continue;
        end
        
        % Get x and y coordinates
        x = boundary(:, 2);
        y = boundary(:, 1);
        
        % Calculate curvature using finite differences
        dx = diff(x);
        dy = diff(y);
        curvature = abs(diff(atan2(dy, dx)));
        
        % Sum curvature
        totalCurvature = totalCurvature + sum(curvature);
        validBoundaryCount = validBoundaryCount + 1;
    end
    
    % Calculate average curvature
    if validBoundaryCount > 0
        avgCurvature = totalCurvature / validBoundaryCount;
        smoothCurvatureVal = 1 - (avgCurvature / pi); % Normalize
    else
        smoothCurvatureVal = 0;
    end
end

function continuousStrokeVal = extractContinuousStroke(skeleton)
    % Remove short branches
    skeleton = bwmorph(skeleton, 'spur', 3);
    
    % Identify endpoints
    endpointsImg = bwmorph(skeleton, 'endpoints');
    endpointCount = sum(endpointsImg(:));
    
    % Measure total skeleton length
    skeletonLength = sum(skeleton(:));
    
    % Define a continuity measure
    if skeletonLength == 0
        continuousStrokeVal = 0;
    else
        continuousStrokeVal = 1 - (endpointCount / skeletonLength);
    end
end

function separateLetters = extractSeparateLetters(bw_noiseRemoved)
    % EXTRACTSEPARATELETTERS  Counts the number of connected components in the
    % binary image (e.g., after noise removal). Each connected component typically
    % corresponds to a separate letter in print handwriting.
    %
    % INPUT:
    %   bw_noiseRemoved : Binary image after noise removal, where letters should
    %                     be clearly separated in print handwriting.
    %
    % OUTPUT:
    %   separateLetters : The total number of connected components (separate letters).

    % Label the connected components in the image
    cc = bwconncomp(bw_noiseRemoved);

    % Return the number of connected components (separate letters)
    separateLetters = cc.NumObjects;
end

function uprightOrientation = extractMinimalSlant(bw_noiseRemoved)
    % EXTRACTMINIMALSLA NT  Measures the slant of the handwriting in the binary
    % image (e.g., after noise removal). A minimal slant would indicate upright
    % print handwriting, close to a vertical orientation.
    %
    % INPUT:
    %   bw_noiseRemoved : Binary image after noise removal, with the print text
    %                     assumed to be near vertical.
    %
    % OUTPUT:
    %   uprightOrientation : The skew angle (in degrees). A value close to 0
    %                        means minimal slant and upright orientation.

    % Find the orientation of the text (skew angle)
    stats = regionprops(bw_noiseRemoved, 'Orientation');

    % If there are multiple components, compute the average orientation
    if numel(stats) > 1
        orientations = [stats.Orientation];
        uprightOrientation = mean(orientations);
    else
        uprightOrientation = stats.Orientation;
    end

    % Return the skew angle, ideally close to 0 for minimal slant
end

function balancedStrokeShapes = extractBalancedStrokeShapes(bw_combined)
    % EXTRACTBALANCEDSTROKESHAPES  Analyzes the stroke shapes in the image,
    % looking for moderate corners and curves. Balanced stroke shapes
    % are typical of print handwriting, with moderate curves and few sharp bends.
    %
    % INPUT:
    %   bw_combined : Binary image after combining horizontal and vertical
    %                 closing, representing the overall structure of printed
    %                 handwriting strokes.
    %
    % OUTPUT:
    %   balancedStrokeShapes : A scalar value that indicates how balanced the
    %                          stroke shapes are, with higher values meaning
    %                          more balanced strokes (moderate curves and corners).

    % Extract the boundary of the shapes
    boundaryImg = bwperim(bw_combined);

    % 1. Find corner points in the boundary
    cornerMask = bwmorph(boundaryImg, 'branchpoints');
    cornerCount = sum(cornerMask(:));

    % 2. Analyze the smoothness of curves (less sharp corners)
    smoothnessMask = bwmorph(boundaryImg, 'endpoints');
    endpointCount = sum(smoothnessMask(:));

    % 3. Calculate a balanced stroke measure
    % We want moderate curves, so we check that there are some corners
    % but not too many sharp ones.
    totalBoundaryPixels = sum(boundaryImg(:));

    if totalBoundaryPixels == 0
        balancedStrokeShapes = 0; % No boundary, no stroke shape
    else
        % Balance stroke shapes based on corner count and smoothness (less endpoints)
        balancedStrokeShapes = 1 - (cornerCount / totalBoundaryPixels) - (endpointCount / totalBoundaryPixels);
    end

end

function linePresence = extractLinePresence(bw_vertical, bw_horizontal)
    % extractLinePresence computes a simple metric indicating the presence
    % of straight lines in the image. It averages the ratio of foreground pixels
    % in the vertical and horizontal images relative to image area.
    %
    % Input:
    %   bw_vertical   : Binary image emphasizing vertical lines.
    %   bw_horizontal : Binary image emphasizing horizontal lines.
    %
    % Output:
    %   linePresence  : A value between 0 and 1 representing line presence.

    totalPixels = numel(bw_vertical);
    verticalScore = sum(bw_vertical(:)) / totalPixels;
    horizontalScore = sum(bw_horizontal(:)) / totalPixels;
    % Average the scores; higher values indicate more straight lines.
    linePresence = (verticalScore + horizontalScore) / 2;
end

function angleCorners = extractRightAngleCorners(bw_combined)
    % extractRightAngleCorners estimates the density of right-angle (or sharp)
    % corners in the image by using MATLAB's built-in corner detector.
    %
    % Input:
    %   bw_combined : Binary image that combines horizontal and vertical closing.
    %
    % Output:
    %   angleCorners: A normalized value reflecting the number of detected corners.
    
    % Use Harris corner detection (or other method) to detect corners.
    detectedCorners = corner(bw_combined, 'Harris');
    % Normalize the corner count by the image area.
    totalPixels = numel(bw_combined);
    angleCorners = size(detectedCorners, 1) / totalPixels;
end

function strokeLengthConsistency = extractStrokeLengthConsistency(bw)
    % extractStrokeLengthConsistency measures how uniform the stroke lengths
    % are in the handwriting. A consistent stroke length suggests regular,
    % uniform block letters.
    %
    % Input:
    %   bw : Binary image of the handwriting (not necessarily skeletonized).
    %
    % Output:
    %   strokeLengthConsistency : A value between 0 and 1. Values closer to 1
    %                             denote more uniform stroke lengths.
    
    % Extract the skeleton to measure stroke lengths.
    skeleton = bwskel(bw);
    cc = bwconncomp(skeleton);
    stats = regionprops(cc, 'Area');
    areas = [stats.Area];
    
    if isempty(areas) || mean(areas) == 0
        strokeLengthConsistency = 0;
    else
        % Compute coefficient of variation as a measure of inconsistency.
        cv = std(areas) / mean(areas);
        % Define consistency as the inverse of variability, normalized to [0,1].
        strokeLengthConsistency = max(0, min(1, 1 - cv));
    end
end

function avgSlantAngle = extractSlantAngle(bw)
    % extractSlantAngle computes the average orientation (in degrees)
    % from connected components in the binary image.
    %
    % Input:
    %   bw - Binary image (typically after noise removal)
    %
    % Output:
    %   avgSlantAngle - Average slant angle (in degrees). Positive values
    %                   indicate counter-clockwise tilt.
    
    % Obtain region properties from connected components
    stats = regionprops(bw, 'Orientation', 'Area');
    if isempty(stats)
        avgSlantAngle = 0;
        return;
    end
    
    % Filter out very small components (noise) by thresholding area
    areas = [stats.Area];
    minAreaThreshold = 0.05 * max(areas);
    validIdx = areas >= minAreaThreshold;
    if ~any(validIdx)
        validIdx = true(1, length(stats));
    end
    
    orientations = [stats(validIdx).Orientation];
    
    % Calculate the average orientation as the slant angle.
    avgSlantAngle = mean(orientations);
end

function tiltUniformity = extractLetterTiltUniformity(bw_inverted)
    % extractLetterTiltUniformity measures how uniform the letter tilt is.
    % It computes the standard deviation of connected component orientations
    % and returns a uniformity score normalized to [0, 1] (1 = perfectly uniform).
    %
    % Input:
    %   bw_inverted - Binary image where text is white on black background.
    %
    % Output:
    %   tiltUniformity - Uniformity score (1 means all letters have almost the same tilt).
    
    % Use regionprops to get orientations (filter with area threshold)
    stats = regionprops(bw_inverted, 'Orientation', 'Area');
    if isempty(stats)
        tiltUniformity = 1;
        return;
    end
    
    areas = [stats.Area];
    minAreaThreshold = 0.05 * max(areas);
    validIdx = areas >= minAreaThreshold;
    if ~any(validIdx)
        validIdx = true(1, length(stats));
    end
    
    orientations = [stats(validIdx).Orientation];
    stdOrient = std(orientations);
    
    % Assume a maximum expected standard deviation of 45°.
    % Uniformity is higher when std is lower.
    tiltUniformity = max(0, 1 - (stdOrient / 45));
end

function verticalStrokeCount = extractVerticalStrokeCount(bw_vertical)
    % extractVerticalStrokeCount counts the number of vertical strokes
    % in the binary image (bw_vertical) that emphasize vertical features.
    % In slanted handwriting, fewer distinct vertical strokes are expected.
    %
    % Input:
    %   bw_vertical - Binary image emphasizing vertical strokes (from morphological operations).
    %
    % Output:
    %   verticalStrokeCount - Total count of connected components that behave as vertical strokes.
    
    cc = bwconncomp(bw_vertical);
    stats = regionprops(cc, 'BoundingBox');
    
    count = 0;
    for i = 1:length(stats)
        bbox = stats(i).BoundingBox; % Format: [x, y, width, height]
        % Consider a component vertical if its height is at least 3 times its width.
        if bbox(4) / bbox(3) > 3
            count = count + 1;
        end
    end
    verticalStrokeCount = count;
end

function cornerCount = extractCornerCount(bw_combined)
    % extractCornerCount detects and counts the number of sharp corners
    % in handwriting using Harris corner detection.
    %
    % Input:
    %   bw_combined - Binary image after noise removal and morphological operations.
    %
    % Output:
    %   cornerCount - Total number of detected corners in the handwriting.
    
    % Apply Harris corner detection with appropriate sensitivity
    corners = corner(bw_combined, 'Harris', 'Sensitivity', 0.15);
    
    % Get the total number of detected corners
    cornerCount = size(corners, 1);
    
    % Normalize by image size if needed
    % totalPixels = numel(bw_combined);
    % cornerCount = cornerCount / sqrt(totalPixels) * 100;
end

function circularity = extractCircularity(bw)
    % extractCircularity calculates the average circularity of letters
    % in the handwriting. Angular handwriting has lower circularity scores.
    %
    % Input:
    %   bw - Binary image of the handwriting.
    %
    % Output:
    %   circularity - Average circularity measure [0-1]. Lower values
    %                 indicate more angular shapes.
    
    % Label connected components
    cc = bwconncomp(bw);
    
    % If no components, return 0 (no circularity)
    if cc.NumObjects == 0
        circularity = 0;
        return;
    end
    
    % Get region properties needed for circularity calculation
    stats = regionprops(cc, 'Area', 'Perimeter');
    
    % Filter out very small components (likely noise)
    areas = [stats.Area];
    perimeters = [stats.Perimeter];
    
    % Remove components that are too small
    minSize = 10; % Minimum area to consider
    validIdx = areas > minSize & perimeters > 0;
    
    if ~any(validIdx)
        circularity = 0;
        return;
    end
    
    areas = areas(validIdx);
    perimeters = perimeters(validIdx);
    
    % Calculate circularity for each component
    % Circularity = 4π × Area / Perimeter²
    % Perfect circle has circularity of 1, more angular shapes are closer to 0
    individualCircularities = 4 * pi * areas ./ (perimeters.^2);
    
    % Some values might be > 1 due to discretization errors, cap at 1
    individualCircularities = min(individualCircularities, 1);
    
    % Calculate the average circularity, weighted by area
    circularity = sum(individualCircularities .* areas) / sum(areas);
end

function cornerDensity = extractCornerDensity(bw_combined)
    % extractCornerDensity calculates the density of corners in handwriting.
    % Angular handwriting has a higher corner density.
    %
    % Input:
    %   bw_combined - Binary image after noise removal and morphological operations.
    %
    % Output:
    %   cornerDensity - A measure of corners per connected component.
    
    % Get connected components
    cc = bwconncomp(bw_combined);
    
    if cc.NumObjects == 0
        cornerDensity = 0;
        return;
    end
    
    % Detect corners using Harris corner detector
    corners = corner(bw_combined, 'Harris', 'Sensitivity', 0.15);
    totalCorners = size(corners, 1);
    
    % Calculate total area of all components
    stats = regionprops(cc, 'Area');
    totalArea = sum([stats.Area]);
    
    if totalArea == 0
        cornerDensity = 0;
        return;
    end
    
    % Calculate corners per unit area
    cornerDensity = totalCorners / sqrt(totalArea);
    
    % Normalize to a reasonable range (0-1)
    cornerDensity = min(cornerDensity / 0.1, 1);
end

function strokeWidthVariation = extractStrokeWidthVariation(bw)
    % extractStrokeWidthVariation detects the variation in stroke width 
    % characteristic of modern calligraphy with thick-thin contrasts.
    %
    % Input:
    %   bw - Binary image of the handwriting
    %
    % Output:
    %   strokeWidthVariation - Measure of stroke width variation [0-1]
    %                         Higher values indicate more thick-thin contrast
    
    % Calculate distance transform for stroke width estimation
    distMap = bwdist(~bw);
    
    % Get non-zero values which represent stroke width at each point
    strokeWidths = distMap(distMap > 0);
    
    if isempty(strokeWidths)
        strokeWidthVariation = 0;
        return;
    end
    
    % Calculate coefficient of variation (CV = standard deviation / mean)
    meanWidth = mean(strokeWidths);
    stdWidth = std(strokeWidths);
    
    if meanWidth == 0
        strokeWidthVariation = 0;
    else
        % Normalize CV to a 0-1 scale (typical calligraphy has high variation)
        cv = stdWidth / meanWidth;
        
        % Map CV to a 0-1 scale, assuming CV values typically range from 0 to 1
        % for handwriting (higher for calligraphy)
        strokeWidthVariation = min(cv / 0.8, 1);
    end
end

function flourishes = extractFlourishes(bw_disk)
    % extractFlourishes detects decorative elements like loops and extended strokes 
    % characteristic of calligraphy flourishes.
    %
    % Input:
    %   bw_disk - Binary image processed with disk-based operations
    %
    % Output:
    %   flourishes - Measure of flourish presence [0-1]
    
    % 1. Skeletonize the image to analyze the structure
    skel = bwskel(bw_disk);
    
    % 2. Find branch points (junctions) and endpoints
    branchPoints = bwmorph(skel, 'branchpoints');
    endPoints = bwmorph(skel, 'endpoints');
    
    numBranchPoints = sum(branchPoints(:));
    numEndPoints = sum(endPoints(:));
    
    % 3. Calculate total skeleton length
    skelLength = sum(skel(:));
    
    if skelLength == 0
        flourishes = 0;
        return;
    end
    
    % 4. Analyze connected components
    cc = bwconncomp(bw_disk);
    stats = regionprops(cc, 'BoundingBox', 'Area', 'Perimeter');
    
    % Calculate expected baseline height
    heights = zeros(length(stats), 1);
    for i = 1:length(stats)
        heights(i) = stats(i).BoundingBox(4);
    end
    medHeight = median(heights(heights > 5)); % Filter out noise
    
    if isempty(medHeight) || medHeight == 0
        flourishes = 0;
        return;
    end
    
    % 5. Identify potential flourishes by:
    %    - Components with high perimeter-to-area ratio (extended strokes)
    %    - Components much larger than median height (extended decorative elements)
    flourishScore = 0;
    for i = 1:length(stats)
        area = stats(i).Area;
        perim = stats(i).Perimeter;
        height = stats(i).BoundingBox(4);
        
        % Skip very small components
        if area < 20
            continue;
        end
        
        % Calculate elongation factor
        elongation = perim^2 / (4 * pi * area);
        
        % Check for height exceeding typical character height
        heightRatio = height / medHeight;
        
        % Add to flourish score if elongated or much taller
        if elongation > 3 || heightRatio > 1.5
            flourishScore = flourishScore + area * (elongation/10 + heightRatio/2);
        end
    end
    
    % 6. Normalize flourish score based on total image area
    totalArea = numel(bw_disk);
    flourishes = min(flourishScore / (totalArea * 0.05), 1);
    
    % 7. Enhance with branch point information
    % More branch points relative to endpoints can indicate decorative elements
    branchToEndRatio = numBranchPoints / max(1, numEndPoints);
    flourishes = flourishes * 0.7 + min(branchToEndRatio/3, 1) * 0.3;
    
    flourishes = min(flourishes, 1);
end

function smoothCurves = extractSmoothCurves(bw_disk)
    % extractSmoothCurves analyzes the smoothness of curves in the handwriting,
    % characteristic of modern calligraphy's flowing style.
    %
    % Input:
    %   bw_disk - Binary image processed with disk-based operations
    %
    % Output:
    %   smoothCurves - Measure of curve smoothness [0-1]
    %                  Higher values indicate smoother curves
    
    % 1. Get boundaries of all objects
    [boundaries, ~] = bwboundaries(bw_disk);
    
    if isempty(boundaries)
        smoothCurves = 0;
        return;
    end
    
    totalCurvature = 0;
    totalPoints = 0;
    
    % 2. Process each boundary to calculate curvature
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        
        % Skip very small boundaries
        if size(boundary, 1) < 10
            continue;
        end
        
        % Get x and y coordinates
        x = boundary(:, 2);
        y = boundary(:, 1);
        
        % Smooth the boundary slightly to reduce pixelation effects
        x = smooth(x, 5);
        y = smooth(y, 5);
        
        % Calculate local curvature using finite differences
        n = length(x);
        curvature = zeros(n, 1);
        
        for i = 2:(n-1)
            % Finite difference approximation of curvature
            dx1 = x(i) - x(i-1);
            dy1 = y(i) - y(i-1);
            dx2 = x(i+1) - x(i);
            dy2 = y(i+1) - y(i);
            
            % Angle change - less change means smoother curve
            angle1 = atan2(dy1, dx1);
            angle2 = atan2(dy2, dx2);
            angleDiff = abs(angle2 - angle1);
            
            % Normalize to [0, π]
            angleDiff = mod(angleDiff + pi, 2*pi) - pi;
            
            curvature(i) = abs(angleDiff);
        end
        
        % Fill in boundary points
        curvature(1) = curvature(2);
        curvature(n) = curvature(n-1);
        
        % Add to total
        totalCurvature = totalCurvature + sum(curvature);
        totalPoints = totalPoints + n;
    end
    
    if totalPoints == 0
        smoothCurves = 0;
    else
        % Calculate average curvature
        avgCurvature = totalCurvature / totalPoints;
        
        % Convert average curvature to smoothness measure
        % Lower curvature means smoother curves
        smoothCurves = max(0, 1 - (avgCurvature / (pi/2)));
    end
end
