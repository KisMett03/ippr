% batchProcessHandwritingImages();

function batchProcessHandwritingImages()

    folderPath = 'Test data/Sloppy';
    
    % Check for valid folder
    if ~exist(folderPath, 'dir')
        error('Folder not found: %s', folderPath);
    end
    
    % Get all image files in the folder
    filePattern = fullfile(folderPath, '*.png');
    pngFiles = dir(filePattern);
    
    filePattern = fullfile(folderPath, '*.jpg');
    jpgFiles = dir(filePattern);
    
    filePattern = fullfile(folderPath, '*.jpeg');
    jpegFiles = dir(filePattern);
    
    filePattern = fullfile(folderPath, '*.bmp');
    bmpFiles = dir(filePattern);
    
    imageFiles = [pngFiles; jpgFiles; jpegFiles; bmpFiles];
    
    if isempty(imageFiles)
        error('No image files found in the folder');
    end
    
    % Initialize feature table
    numImages = length(imageFiles);
    
    % Column names for features
    featureNames = {'Filename', 'AvgHeight', 'AvgWidth', 'StrokeWidth', ...
                   'BaselineAngle', 'LineSpacing', 'WordSpacing', 'Connectivity', ...
                   'HeightConsistency', 'TextDensity', 'SlantAngle', ...
                   'IsAligned', 'StrokeConsistency', 'IrregularLetterSize', ...
                   'SmoothCurves', 'PrimaryStyle', 'PrimaryConfidence'};
    
    % Add style label if provided
    if nargin > 2 && ~isempty(styleLabel)
        featureNames = [featureNames, 'ActualStyle'];
    end
    
    % Initialize data arrays
    featureData = cell(numImages, length(featureNames));
    
    % Process each image
    fprintf('Processing %d images from %s\n', numImages, folderPath);
    
    for i = 1:numImages
        try
            % Get the full filename
            filename = fullfile(folderPath, imageFiles(i).name);
            fprintf('Processing image %d/%d: %s\n', i, numImages, imageFiles(i).name);
            
            % Load the image
            img = imread(filename);
            
            % Convert to grayscale if it's RGB
            if size(img, 3) == 3
                img_gray = rgb2gray(img);
            else
                img_gray = img;
            end
            
            % Enhance contrast using CLAHE for better feature detection
            img_enhanced = adapthisteq(img_gray, 'ClipLimit', 0.02);
            
            % Apply sharpening to enhance handwriting strokes
            sharpened = imsharpen(img_enhanced, 'Radius', 1, 'Amount', 0.5);
            
            % Step 2: Binarize using adaptive thresholding
            bw = imbinarize(sharpened, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.4);
            
            % Step 3: Invert the binary image
            bw = ~bw;
            
            % Step 4: Remove small objects (noise)
            bw = bwareaopen(bw, 20);
            
            % Original binary image (save for specific features)
            bw_original = bw;
            
            % Step 5: Apply disk-based opening for curved strokes preservation
            se_disk = strel('disk', 1);
            bw_disk = imopen(bw, se_disk);
            
            % Step 6: Apply horizontal line-based closing for baseline features
            se_horizontal = strel('line', 5, 0);
            bw_horizontal = imclose(bw, se_horizontal);
            
            % Step 7: Apply vertical line-based closing for letter connections
            se_vertical = strel('line', 3, 90);
            bw_vertical = imclose(bw, se_vertical);
            
            % Step 8: Combined processing for skeleton extraction
            bw_combined = bw;
            bw_combined = imclose(bw_combined, strel('line', 3, 0));
            bw_combined = imclose(bw_combined, strel('line', 3, 90));
            skeleton = bwskel(bw_combined);
            
            % Extract features
            [avgHeight, avgWidth] = detectLetterSize(skeleton);
            avgStrokeWidth = detectStrokeWidth(bw_original);
            baselineAngle = detectBaselineAngle(bw_horizontal);
            avgLineSpacing = detectLineSpacing(bw_vertical);
            avgWordSpacing = detectWordSpacing(bw_horizontal);
            connectivity = detectLetterConnectivity(bw_combined);
            heightConsistency = detectLetterHeightConsistency(skeleton);
            textDensity = detectTextDensity(bw_original);
            slantAngle = detectSlantAngle(bw_original);
            isAligned = detectBaselineAlignment(bw_horizontal);
            strokeConsistency = detectStrokeConsistency(bw_original);
            irregularLetterSize = detectIrregularLetterSize(skeleton);
            smoothCurves = detectSmoothCurves(bw_disk);
            
            % Classify style
            styles = classifyHandwritingStyle(slantAngle, connectivity, avgWordSpacing, ...
                smoothCurves, heightConsistency, avgHeight, irregularLetterSize, ...
                strokeConsistency, isAligned);
            
            % Store features
            featureData{i, 1} = imageFiles(i).name;
            featureData{i, 2} = avgHeight;
            featureData{i, 3} = avgWidth;
            featureData{i, 4} = avgStrokeWidth;
            featureData{i, 5} = baselineAngle;
            featureData{i, 6} = avgLineSpacing;
            featureData{i, 7} = avgWordSpacing;
            featureData{i, 8} = connectivity;
            featureData{i, 9} = heightConsistency;
            featureData{i, 10} = textDensity;
            featureData{i, 11} = slantAngle;
            featureData{i, 12} = isAligned;
            featureData{i, 13} = strokeConsistency;
            featureData{i, 14} = irregularLetterSize;
            featureData{i, 15} = smoothCurves;
            featureData{i, 16} = styles.primaryStyle;
            featureData{i, 17} = styles.primaryConfidence;
            
            % Add style label if provided
            if nargin > 2 && ~isempty(styleLabel)
                featureData{i, 18} = styleLabel;
            end
            
        catch ex
            warning('Error processing %s: %s', imageFiles(i).name, ex.message);
            % Fill with NaN or empty values for this image
            featureData{i, 1} = imageFiles(i).name;
            featureData{i, 2} = NaN;
            featureData{i, 3} = NaN;
            % Continue filling NaNs for all numeric features
            for j = 4:15
                featureData{i, j} = NaN;
            end
            featureData{i, 16} = 'error';
            featureData{i, 17} = NaN;
            
            if nargin > 2 && ~isempty(styleLabel)
                featureData{i, 18} = styleLabel;
            end
        end
    end
    
    % Convert to table and save to CSV
    featureTable = cell2table(featureData, 'VariableNames', featureNames);
    writetable(featureTable, 'sloppy.csv');
    
    fprintf('Processing complete');
end