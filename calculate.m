% Handwriting Style Feature Analysis
% This script analyzes key features for different handwriting styles

% Load the CSV files
italicData = readtable('italic.csv');
printData = readtable('print.csv');
roundData = readtable('round.csv');
sloppyData = readtable('sloppy.csv');

% 1. ITALIC STYLE - Key features: SlantAngle and AvgHeight
fprintf('=== ITALIC HANDWRITING FEATURES ===\n');
% Calculate statistics for SlantAngle
italicSlant = italicData.SlantAngle;
italicSlant = italicSlant(~isnan(italicSlant)); % Remove NaN values
fprintf('Letter Slant (degrees):\n');
fprintf('  Median: %.2f\n', median(italicSlant));
fprintf('  Q1: %.2f\n', quantile(italicSlant, 0.25));
fprintf('  Q3: %.2f\n', quantile(italicSlant, 0.75));

% Calculate statistics for AvgHeight
italicHeight = italicData.AvgHeight;
italicHeight = italicHeight(~isnan(italicHeight)); % Remove NaN values
fprintf('Letter Height (pixels):\n');
fprintf('  Median: %.2f\n', median(italicHeight));
fprintf('  Q1: %.2f\n', quantile(italicHeight, 0.25));
fprintf('  Q3: %.2f\n', quantile(italicHeight, 0.75));

% 2. PRINT STYLE - Key features: AvgWordSpacing, SmoothCurves, HeightConsistency
fprintf('\n=== PRINT (BLOCK) HANDWRITING FEATURES ===\n');
% Calculate statistics for AvgWordSpacing
printSpacing = printData.WordSpacing;
printSpacing = printSpacing(~isnan(printSpacing));
fprintf('Letter Spacing (pixels):\n');
fprintf('  Median: %.2f\n', median(printSpacing));
fprintf('  Q1: %.2f\n', quantile(printSpacing, 0.25));
fprintf('  Q3: %.2f\n', quantile(printSpacing, 0.75));

% Calculate statistics for SmoothCurves (inverse = stroke sharpness)
printSharpness = 1 - printData.SmoothCurves;
printSharpness = printSharpness(~isnan(printSharpness));
fprintf('Stroke Sharpness (0-1):\n');
fprintf('  Median: %.2f\n', median(printSharpness));
fprintf('  Q1: %.2f\n', quantile(printSharpness, 0.25));
fprintf('  Q3: %.2f\n', quantile(printSharpness, 0.75));

% Calculate statistics for HeightConsistency
printUniformity = printData.HeightConsistency;
printUniformity = printUniformity(~isnan(printUniformity));
fprintf('Character Uniformity (0-1):\n');
fprintf('  Median: %.2f\n', median(printUniformity));
fprintf('  Q1: %.2f\n', quantile(printUniformity, 0.25));
fprintf('  Q3: %.2f\n', quantile(printUniformity, 0.75));

% 3. ROUNDED STYLE - Key features: SmoothCurves and IrregularLetterSize
fprintf('\n=== ROUNDED HANDWRITING FEATURES ===\n');
% Calculate statistics for SmoothCurves
roundCurves = roundData.SmoothCurves;
roundCurves = roundCurves(~isnan(roundCurves));
fprintf('Smooth Curves (0-1):\n');
fprintf('  Median: %.2f\n', median(roundCurves));
fprintf('  Q1: %.2f\n', quantile(roundCurves, 0.25));
fprintf('  Q3: %.2f\n', quantile(roundCurves, 0.75));

% Calculate statistics for Consistent Spacing (inverse of IrregularLetterSize)
roundConsistency = 1 - roundData.IrregularLetterSize;
roundConsistency = roundConsistency(~isnan(roundConsistency));
fprintf('Consistent Spacing (0-1):\n');
fprintf('  Median: %.2f\n', median(roundConsistency));
fprintf('  Q1: %.2f\n', quantile(roundConsistency, 0.25));
fprintf('  Q3: %.2f\n', quantile(roundConsistency, 0.75));

% 4. SLOPPY STYLE - Key features: IrregularLetterSize and (1-HeightConsistency)
fprintf('\n=== SLOPPY HANDWRITING FEATURES ===\n');
% Calculate statistics for IrregularLetterSize
sloppyIrregularity = sloppyData.IrregularLetterSize;
sloppyIrregularity = sloppyIrregularity(~isnan(sloppyIrregularity));
fprintf('Irregular Letter Size (0-1):\n');
fprintf('  Median: %.2f\n', median(sloppyIrregularity));
fprintf('  Q1: %.2f\n', quantile(sloppyIrregularity, 0.25));
fprintf('  Q3: %.2f\n', quantile(sloppyIrregularity, 0.75));

% Calculate statistics for Height Inconsistency (inverse of HeightConsistency)
sloppyInconsistency = 1 - sloppyData.HeightConsistency;
sloppyInconsistency = sloppyInconsistency(~isnan(sloppyInconsistency));
fprintf('Inconsistent Heights (0-1):\n');
fprintf('  Median: %.2f\n', median(sloppyInconsistency));
fprintf('  Q1: %.2f\n', quantile(sloppyInconsistency, 0.25));
fprintf('  Q3: %.2f\n', quantile(sloppyInconsistency, 0.75));