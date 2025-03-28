function main()
    % Create a figure for the GUI
    fig = uifigure('Position', [100, 100, 1000, 700], 'Name', 'Handwriting Analysis');

    % Create a button to select the image
    btnSelectImage = uibutton(fig, 'Position', [50, 650, 150, 30], 'Text', 'Select Image', 'ButtonPushedFcn', @(btn, event) selectImage());

    % Create a button to run the analysis
    btnAnalyze = uibutton(fig, 'Position', [220, 650, 150, 30], 'Text', 'Analyze', 'ButtonPushedFcn', @(btn, event) analyzeImage());

    % Create panels for images
    originalPanel = uipanel(fig, 'Title', 'Original Image', 'Position', [50, 400, 450, 230]);
    processedPanel = uipanel(fig, 'Title', 'Processed Image', 'Position', [500, 400, 450, 230]);

    % Create panel for feature results
    featurePanel = uipanel(fig, 'Title', 'Handwriting Features with nomalise values', 'Position', [50, 150, 450, 230]);

    % Create text area for feature results display
    featureText = uitextarea(featurePanel, 'Position', [10, 10, 430, 190], 'Editable', 'off');
    featureText.FontName = 'Consolas'; % Monospaced font for better alignment
    featureText.FontSize = 16;

    % Create panel for style classification results
    stylePanel = uipanel(fig, 'Title', 'Style Classification', 'Position', [500, 150, 450, 230]);

    % Create text area for style classification display
    styleText = uitextarea(stylePanel, 'Position', [10, 10, 430, 190], 'Editable', 'off');
    styleText.FontName = 'Consolas';
    styleText.FontSize = 16;

    % Create axes for displaying images
    originalAxes = uiaxes(originalPanel);
    originalAxes.Position = [10, 10, 430, 190];
    originalAxes.XTick = [];
    originalAxes.YTick = [];

    processedAxes = uiaxes(processedPanel);
    processedAxes.Position = [10, 10, 430, 190];
    processedAxes.XTick = [];
    processedAxes.YTick = [];

    % Declare the image path for analysis
    imgPath = '';
    inputImg = [];

    % Declare variables to hold the processed images
    processedImages = struct();

    % Function to select image
    function selectImage()
        [filename, pathname] = uigetfile({'*.png;*.jpg;*.jpeg;*.bmp;*.tif;*.tiff', 'Image Files'}, 'Select an image file');

        if filename ~= 0
            imgPath = fullfile(pathname, filename);
            disp(['Selected Image: ', imgPath]);

            % Display the selected image
            inputImg = imread(imgPath);
            imshow(inputImg, 'Parent', originalAxes);
            title(originalAxes, 'Original Image');

            % Clear processed images when new image is selected
            cla(processedAxes);

            % Clear feature and style text areas
            featureText.Value = '';
            styleText.Value = '';
        else
            disp('No file selected');
        end

    end

    % Function to run analysis
    function analyzeImage()

        if isempty(imgPath)
            uialert(fig, 'Please select an image first.', 'Error', 'Icon', 'error');
            return;
        end

        % Load and process the image
        [bw, bw_inverted, bw_noiseRemoved, bw_disk, bw_horizontal, bw_vertical, skeleton] = preprocess(inputImg);

        % Store processed images in the struct
        processedImages.bw = bw;
        processedImages.bw_inverted = bw_inverted;
        processedImages.bw_noiseRemoved = bw_noiseRemoved;
        processedImages.bw_disk = bw_disk;
        processedImages.bw_horizontal = bw_horizontal;
        processedImages.bw_vertical = bw_vertical;
        processedImages.skeleton = skeleton;

        % Save the processed images to files
        savePath = fullfile(fileparts(imgPath), 'processed_images');

        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end

        % Save each processed image with a descriptive filename
        imwrite(bw, fullfile(savePath, 'binary.png'));
        imwrite(bw_inverted, fullfile(savePath, 'inverted.png'));
        imwrite(bw_noiseRemoved, fullfile(savePath, 'noise_removed.png'));
        imwrite(bw_disk, fullfile(savePath, 'morphology_disk.png'));
        imwrite(bw_horizontal, fullfile(savePath, 'horizontal_elements.png'));
        imwrite(bw_vertical, fullfile(savePath, 'vertical_elements.png'));
        imwrite(skeleton, fullfile(savePath, 'skeleton.png'));

        % Display confirmation message
        disp('Processed images saved to: ' + string(savePath));
        % Display processed images
        imshow(bw_noiseRemoved, 'Parent', processedAxes);
        title(processedAxes, 'Processed Image');

        % Call feature extraction function
        disp('Extracting features...');

        % Create bw_combined for feature extraction
        bw_combined = bw_noiseRemoved;

        % Call the feature extraction function
        [group1Features, group2Features, group3Features] = feature_ext(bw, bw_inverted, bw_noiseRemoved, bw_disk, bw_horizontal, bw_vertical, skeleton, bw_combined);

        % Display extracted features in the text area
        displayFeatures(group1Features, group2Features, group3Features);

        % Perform style classification based on extracted features
        classifyHandwritingStyle(group1Features, group2Features, group3Features);
    end

    % Function to display extracted features
    function displayFeatures(group1, group2, group3)
        % Create a string cell array to hold feature results
        featLines = cell(30, 1);
        lineIdx = 1;

        % Group 1: Cursive & Print Handwriting
        featLines{lineIdx} = '=== CURSIVE FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Connected Components: %.2f', group1.cursive.connectedComponents); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Smooth Curvature: %.2f', group1.cursive.smoothCurvature); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Continuous Stroke: %.2f', group1.cursive.continuousStroke); lineIdx = lineIdx + 1;
        featLines{lineIdx} = ' '; lineIdx = lineIdx + 1;

        % Group 1: Print Features
        featLines{lineIdx} = '=== PRINT FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Separate Letters: %.2f', group1.print.separateLetters); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Upright Orientation: %.2f', group1.print.uprightOrientation); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Balanced Stroke Shapes: %.2f', group1.print.balancedStrokeShapes); lineIdx = lineIdx + 1;
        featLines{lineIdx} = ' '; lineIdx = lineIdx + 1;

        % Group 2: Block Letters & Slanted Handwriting
        featLines{lineIdx} = '=== BLOCK LETTER FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Line Presence: %.2f', group2.blockLetters.linePresence); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Stroke Length Consistency: %.2f', group2.blockLetters.strokeLengthConsistency); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Dominant HV Lines: %.2f', group2.blockLetters.hvDominance); lineIdx = lineIdx + 1;
        featLines{lineIdx} = ' '; lineIdx = lineIdx + 1;

        % Group 2: Slanted Features
        featLines{lineIdx} = '=== SLANTED FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Avg Slant Angle: %.2f ', group2.slanted.avgSlantAngle); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Letter Tilt Uniformity: %.2f', group2.slanted.letterTiltUniformity); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Vertical Stroke Count: %.2f', group2.slanted.verticalStrokeCount); lineIdx = lineIdx + 1;
        featLines{lineIdx} = ' '; lineIdx = lineIdx + 1;

        % Group 3: Angular Features
        featLines{lineIdx} = '=== ANGULAR FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Edge Orientation Variance: %.2f', group3.angular.edgeOrientationVariance); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Circularity: %.2f', group3.angular.circularity); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Corner Density: %.2f', group3.angular.cornerDensity); lineIdx = lineIdx + 1;
        featLines{lineIdx} = ' '; lineIdx = lineIdx + 1;

        % Group 3: Calligraphy Features
        featLines{lineIdx} = '=== CALLIGRAPHY FEATURES ==='; lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Stroke Width Variation: %.2f', group3.calligraphy.strokeWidthVariation); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Flourishes: %.2f', group3.calligraphy.flourishes); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Smooth Curves: %.2f', group3.calligraphy.smoothCurves); lineIdx = lineIdx + 1;
        featLines{lineIdx} = sprintf('Letter Extension: %.2f', group3.calligraphy.extendedAscDescMeasure); lineIdx = lineIdx + 1;

        % Remove any empty cells (which would be []) before assignment
        featLines = featLines(~cellfun(@isempty, featLines));
        featureText.Value = featLines;
    end

    % Function to classify handwriting style based on features
    function classifyHandwritingStyle(group1, group2, group3)
        % Compute style scores based on normalized features

        styleScores = struct();

        % Cursive score: higher connectivity (closer to 1), smooth curvature, and continuous stroke indicate cursive.
        styleScores.cursive = group1.cursive.connectedComponents * 0.4 + ...
            group1.cursive.smoothCurvature * 0.3 + ...
            group1.cursive.continuousStroke * 0.3;

        % Print score: higher separate letters, upright orientation, and balanced stroke shapes indicate print.
        styleScores.print = group1.print.separateLetters * 0.3 + ...
            group1.print.uprightOrientation * 0.3 + ...
            group1.print.balancedStrokeShapes * 0.4;

        % Block letters score: marked by strong line presence, consistent stroke length, and dominant horizontal/vertical lines.
        styleScores.block = group2.blockLetters.linePresence * 0.35 + ...
            group2.blockLetters.strokeLengthConsistency * 0.4 + ...
            group2.blockLetters.hvDominance * 0.25;

        % Slanted score: a larger avg slant angle (normalized so that 1 means maximum slant) plus uniform tilt
        % and fewer vertical strokes (thus we use 1 - verticalStrokeCount).
        styleScores.slanted = group2.slanted.avgSlantAngle * 0.5 + ...
            group2.slanted.letterTiltUniformity * 0.3 + ...
            (1 - group2.slanted.verticalStrokeCount) * 0.2;

        % Angular score: high edge orientation variance, low circularity (hence 1-circularity),
        % and high corner density indicate an angular style.
        styleScores.angular = group3.angular.edgeOrientationVariance * 0.3 + ...
            (1 - group3.angular.circularity) * 0.4 + ...
            group3.angular.cornerDensity * 0.25;

        % Calligraphy score: high stroke width variation, decorative flourishes,
        % smooth curves, and extended ascenders/descenders.
        styleScores.calligraphy = group3.calligraphy.strokeWidthVariation * 0.3 + ...
            group3.calligraphy.flourishes * 0.25 + ...
            group3.calligraphy.smoothCurves * 0.25 + ...
            group3.calligraphy.extendedAscDescMeasure * 0.2;

        % Gather style names and scores for comparison
        styleNames = fieldnames(styleScores);
        styleValues = zeros(length(styleNames), 1);

        for i = 1:length(styleNames)
            styleValues(i) = styleScores.(styleNames{i});
        end

        % Sort styles by descending score
        [sortedScores, sortIndices] = sort(styleValues, 'descend');

        % Build output lines for style display
        styleLines = cell(7, 1);
        styleLines{1} = 'HANDWRITING STYLE ANALYSIS:';
        styleLines{2} = '========================';

        for i = 1:min(3, length(styleNames))
            styleIdx = sortIndices(i);
            styleLines{i + 2} = sprintf('%d. %s (%.2f%%)', i, upper(styleNames{styleIdx}), sortedScores(i) * 100);
        end

        % Add primary characteristic description based on the highest-scoring style
        primaryStyle = styleNames{sortIndices(1)};

        switch primaryStyle
            case 'cursive'
                styleLines{6} = 'Primary characteristics: Connected letters, flowing strokes,';
                styleLines{7} = 'continuous writing with minimal pen lifts.';
            case 'print'
                styleLines{6} = 'Primary characteristics: Separated letters, upright orientation,';
                styleLines{7} = 'balanced stroke shapes similar to printed text.';
            case 'block'
                styleLines{6} = 'Primary characteristics: Strong straight lines, right-angle corners,';
                styleLines{7} = 'uniform stroke length resembling block capitals.';
            case 'slanted'
                styleLines{6} = 'Primary characteristics: Consistent letter tilt, fewer vertical strokes,';
                styleLines{7} = 'uniform angle across characters.';
            case 'angular'
                styleLines{6} = 'Primary characteristics: Sharp corners, low circularity,';
                styleLines{7} = 'high density of angular features in writing.';
            case 'calligraphy'
                styleLines{6} = 'Primary characteristics: Thick-thin stroke variation, decorative elements,';
                styleLines{7} = 'smooth curves and flowing strokes.';
        end

        % Update the style text area with the results
        styleText.Value = styleLines;
    end

end
