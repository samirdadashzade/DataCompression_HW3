[I,map]=imread('image.gif');
G=ind2gray(I,map);
% imagesc(I); colormap(map);
% imagesc(G); colormap(gray); 

% a) Compute the entropy of the image.
E = entropy(G);
G = double(G);

% b) Uniform
dFirst = min(min(G));
dLast = max(max(G));

[dLevels, rLevels] = calcUniformIntervals(dFirst, dLast, 8);
[GuPrime, GuHat] = quanDequantMatrix(G, dLevels, rLevels);
E_Gu = entropy(uint8(GuPrime));
snrUniform = snr(G, G - GuHat);
% imagesc(GuHat); colormap(gray); 

% c) Semi-Unifrom
[dLevelsSu, rLevelsSu] = calcSemiUniformIntervals(G, dFirst, dLast, 8);
[GsuPrime, GsuHat] = quanDequantMatrix(G, dLevelsSu, rLevelsSu);
E_Gsu = entropy(uint8(GsuPrime));
snrSemiUniform = snr(G, G - GsuHat);
% imagesc(GsuHat); colormap(gray); 

% d) Max-Lloyd
[dLevelsOp, rLevelsOp] = calcOptimalIntervals(G, dFirst, dLast, 8);
[GopPrime, GopHat] = quanDequantMatrix(G, dLevelsOp, rLevelsOp);
E_Gop = entropy(uint8(GopPrime));
snrMaxLloyd = snr(G, G - GopHat);
% imagesc(GopHat); colormap(gray); 

% e) Comparison
label = ["Uniform"; "Semi-Uniform"; "Max-Lloyd"];
entropies = [E_Gu; E_Gsu; E_Gop];
snr = [snrUniform; snrSemiUniform; snrMaxLloyd];
comparisonTable = table(label, entropies, snr);

% imagesc(GuHat); colormap(gray);
% imagesc(GsuHat); colormap(gray);
% imagesc(GopHat); colormap(gray);

function [dLevels, rLevels] = calcUniformIntervals(dFirst, dLast, levelSize)
    dLevels = zeros(1, levelSize + 1, 'double');
    rLevels = zeros(1, levelSize, 'double');
    dLevels(1) = dFirst;
    dLevels(levelSize + 1) = dLast;
    
    delta = (dLast - dFirst)/levelSize;
    for n = 2:levelSize
        dLevels(n) = dLevels(n-1) + delta;
    end
    
    for n = 1:levelSize
        rLevels(n) = (dLevels(n) + dLevels(n+1))/2;
    end
end

function [dLevels, rLevels] = calcSemiUniformIntervals(matrix, dFirst, dLast, levelSize)
    dLevels = zeros(1, levelSize + 1, 'double');
    rLevels = zeros(1, levelSize, 'double');
    dLevels(1) = dFirst;
    dLevels(levelSize + 1) = dLast;
    
    delta = (dLast - dFirst)/levelSize;
    for n = 2:levelSize
        dLevels(n) = dLevels(n-1) + delta;
    end
    
    for n = 1:levelSize
        rLevels(n) = averageBetween(matrix, dLevels(n), dLevels(n+1));
    end
end

function [dLevels, rLevels] = calcOptimalIntervals(matrix, dFirst, dLast, levelSize)
    dLevels = zeros(1, levelSize + 1);
    rLevels = zeros(1, levelSize);
    dLevels(1) = dFirst;
    dLevels(levelSize + 1) = dLast;
    
    delta = (dLast - dFirst)/levelSize;
    for n = 2:levelSize
        dLevels(n) = dLevels(n-1) + delta;
    end

    dLevelsNew = dLevels;
    isFirst = true;

    while isFirst || ~hasConverged(dLevels, dLevelsNew)
        isFirst = false;
        dLevels = dLevelsNew;
        
        for n = 1:levelSize
            rLevels(n) = averageBetween(matrix, dLevels(n), dLevels(n+1));
        end

        for n = 2:levelSize
            dLevelsNew(n) = (rLevels(n-1) + rLevels(n)) / 2;
        end
    end
    
    
end

function converged = hasConverged(arrayA, arrayB)
    diff = 0;

    for i = 1:length(arrayA)
        diff = diff + abs(arrayA(i) - arrayB(i));
    end

    converged = diff == 0;
end

function average = averageBetween(matrix, minValue, maxValue)
    matrixSize = size(matrix);
    sum = double(0);
    count = 0;

    for i = 1:matrixSize(1)
        for j = 1:matrixSize(2)
            if matrix(i, j) >= minValue && matrix(i, j) < maxValue || matrix(i, j) == maxValue
                sum = sum + double(matrix(i, j));
                count = count + 1;
            end
        end
    end   
    average = sum/count;
end

function [quantized, dequantized] = quanDequantMatrix(matrix, dLevels, rLevels)
    matrixSize = size(matrix);
    quantized = zeros(matrixSize(1), matrixSize(2), 'double');
    dequantized = zeros(matrixSize(1), matrixSize(2), 'double');

    for i = 1:matrixSize(1)
        for j = 1:matrixSize(2)
            dIndex = findLevelIndex(dLevels, matrix(i, j));
            quantized(i, j) = dIndex;
            dequantized(i, j) = rLevels(dIndex);
        end
    end    
end

function index = findLevelIndex(levels, value)
    for i = 1:length(levels)-1
        if value == levels(length(levels))
            index = length(levels) - 1;
            break;
        elseif value >= levels(i) && value < levels(i+1)
            index = i;
            break;
        end
    end
end
