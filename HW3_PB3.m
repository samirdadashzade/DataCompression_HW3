[I,map]=imread('image.gif');
G=ind2gray(I,map);
% imagesc(I); colormap(map);
% imagesc(G); colormap(gray); 

X = [0, 0.5, 0.58, 0.69, 0.90, 0.94, 0.97, 1.03, 1.04, 1.27, 1.3, 1.38, 1.4, 1.43, 1.6, 1.7, 1.96, 1.97, 1.99, 2.22 2.27, 2.42, 2.42, 2.57, 2.61, 2.63, 2.83, 3, 3.49, 3.54, 3.66, 3.99];

% a) 
rLevels = calcKmeanIntervals(X, 8);
[quantized, dequantized] = quanDequantArray(X, rLevels);

resultTable = table(X', quantized', dequantized');
resultTable.Properties.VariableNames = ["X"; "IX"; "X_hat"];

mse = immse(X, dequantized);

% b)
G = double(G);
matrixSize = size(G);
imageArray = G(:);
rLevelsMatrix = calcKmeanIntervals(imageArray, 8);
[quantized, dequantized] = quanDequantArray(imageArray, rLevelsMatrix);
quantizedMatrix = reshape(quantized, matrixSize(1), matrixSize(2));
dequantizedMatrix = reshape(dequantized, matrixSize(1), matrixSize(2));

snrKmeans = snr(G, G - dequantizedMatrix);
% imagesc(dequantizedMatrix); colormap(gray);


function rLevels = calcKmeanIntervals(array, levelSize)
    len = length(array);
    clusters = zeros(1, len, 'double');
    rLevels = zeros(1, levelSize, 'double');
    delta = floor(len/levelSize);

    sortedArray = sort(array);
    for i = 1:levelSize
        rLevels(i) = sortedArray(1 + delta * (i - 1));
    end
        
    rLevelsNew = rLevels;
    isFirst = true;

    while isFirst || ~hasConverged(rLevels, rLevelsNew)
        isFirst = false;
        rLevels = rLevelsNew;
        
        disp(rLevels);

        for i = 1:len
            clusters(i) = findMinCluster(array(i), rLevels);
        end
        
        for i = 1:length(rLevels)
            rLevelsNew(i) = calcClusterMean(array, clusters, i);
        end
    end
end

function minIndex = findMinCluster(value, rLevels)
    minIndex = -1;
    minValue = intmax("uint8");

    for i = 1:length(rLevels)
        if abs(value - rLevels(i)) < minValue
            minIndex = i;
            minValue = abs(value - rLevels(i));
        end
    end
end

function clusterMean = calcClusterMean(array, clusters, clusterIndex)
    sum = double(0);
    count = 0;

    for i = 1:length(array)
        if clusters(i) == clusterIndex
            sum = sum + array(i);
            count = count + 1;
        end
    end
    
    if count == 0
        clusterMean = 0;
    else 
        clusterMean = sum/count;
    end
end

function converged = hasConverged(arrayA, arrayB)
    diff = 0;

    for i = 1:length(arrayA)
        diff = diff + abs(arrayA(i) - arrayB(i));
    end

    converged = diff == 0;
end

function [quantized, dequantized] = quanDequantArray(array, rLevels)
    arraySize = length(array);
    quantized = zeros(1, arraySize, 'double');
    dequantized = zeros(1, arraySize, 'double');

    for i = 1:arraySize
        rIndex = findMinCluster(array(i), rLevels);
        quantized(i) = rIndex;
        dequantized(i) = rLevels(rIndex);
    end    
end
