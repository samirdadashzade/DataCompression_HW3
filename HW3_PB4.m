[I,map]=imread('image.gif');
G=ind2gray(I,map);
% imagesc(I); colormap(map);
% imagesc(G); colormap(gray); 

blockSize = 8;

G = double(G);
dG = blockproc(G, [blockSize blockSize], @(blkStruct) dct2(blkStruct.data));
[dcQuantized, dcDequantized] = quantizeDC(dG, blockSize);

matrixSize = size(G);
% imageTotalPixels = width * height
imageTotalPixels = matrixSize(1) * matrixSize(2);
numberOfTotalBlocks = matrixSize(1) * matrixSize(2) / blockSize / blockSize;
% dcTotalBits = numberOfTotalBlocks * bits/dc
dcTotalBits =  numberOfTotalBlocks * log2(8);

% b)
diagonals10 = [2; 3; 4];
[acQuantized10, acDequantized10] = quantizeAC(dG, blockSize, diagonals10);
dequantizedImage10 = dequantizeImage(dG, blockSize, dcDequantized, acDequantized10, diagonals10);
G10 = blockproc(dequantizedImage10, [blockSize blockSize], @(blkStruct) idct2(blkStruct.data));

snrG10 = snr(G, G - G10);

% acTotalBits10 = numberOfTotalBlocks * numberOfDiagonalTerms * bits/ac
acTotalBits10 = numberOfTotalBlocks * 9 * log2(4);
quantizedTotalBits10 = dcTotalBits + acTotalBits10;
bitrate10 = quantizedTotalBits10 / imageTotalPixels; % bits/px
compressionRatio10 = imageTotalPixels * 8 / quantizedTotalBits10;
% imagesc(G10); colormap(gray);

% c)
diagonals6 = [2; 3;];
[acQuantized6, acDequantized6] = quantizeAC(dG, blockSize, diagonals6);
dequantizedImage6 = dequantizeImage(dG, blockSize, dcDequantized, acDequantized6, diagonals6);
G6 = blockproc(dequantizedImage6, [blockSize blockSize], @(blkStruct) idct2(blkStruct.data));

snrG6 = snr(G, G - G6);

% acTotalBits6 = numberOfTotalBlocks * numberOfDiagonalTerms * bits/ac
acTotalBits6 = numberOfTotalBlocks * 5 * log2(4);
quantizedTotalBits6 = dcTotalBits + acTotalBits6;
bitrate6 = quantizedTotalBits6 / imageTotalPixels; % bits/px
compressionRatio6 = imageTotalPixels * 8 / quantizedTotalBits6;
% imagesc(G6); colormap(gray);

% d)
diagonals3 = [2;];
[acQuantized3, acDequantized3] = quantizeAC(dG, blockSize, diagonals3);
dequantizedImage3 = dequantizeImage(dG, blockSize, dcDequantized, acDequantized3, diagonals3);
G3 = blockproc(dequantizedImage3, [blockSize blockSize], @(blkStruct) idct2(blkStruct.data));

snrG3 = snr(G, G - G3);

% acTotalBits3 = numberOfTotalBlocks * numberOfDiagonalTerms * bits/ac
acTotalBits3 = numberOfTotalBlocks * 2 * log2(4);
quantizedTotalBits3 = dcTotalBits + acTotalBits3;
bitrate3 = quantizedTotalBits3 / imageTotalPixels; % bits/px
compressionRatio3 = imageTotalPixels * 8 / quantizedTotalBits3;
% imagesc(G3); colormap(gray);

% e)
dequantizedImage1 = dequantizeImage(dG, blockSize, dcDequantized, [], []);
G1= blockproc(dequantizedImage1, [blockSize blockSize], @(blkStruct) idct2(blkStruct.data));

snrG1 = snr(G, G - G1);

% acTotalBits1 = numberOfTotalBlocks * numberOfDiagonalTerms * bits/ac
quantizedTotalBits1 = dcTotalBits;
bitrate1 = quantizedTotalBits1 / imageTotalPixels; % bits/px
compressionRatio1 = imageTotalPixels * 8 / quantizedTotalBits1;
% imagesc(G1); colormap(gray);


function [dcQuantized, dcDequantized] = quantizeDC(matrix, blockSize)
    matrixSize = size(matrix);
    level1BlockSize = matrixSize(1)/blockSize;
    level2BlockSize = matrixSize(2)/blockSize;
    dcTermsSize = level1BlockSize * level2BlockSize; 
    dcTerms = zeros(1, dcTermsSize, 'double');
    
    for i = 1:level1BlockSize
        for j = 1:level2BlockSize
            value = matrix((i - 1) * blockSize + 1, (j - 1) * blockSize + 1);
            key = (i - 1) * level2BlockSize  + j;
            dcTerms(key) = value;
        end
    end

    dFirst = floor(min(dcTerms));
    dLast = ceil(max(dcTerms) + power(10, -6));
    [dLevels, rLevels] = calcUniformIntervals(dFirst, dLast, 8);
    [dcQuantized, dcDequantized] = quanDequantArray(dcTerms, dLevels, rLevels);
end

function [acQuantized, acDequantized] = quantizeAC(matrix, blockSize, diagonals)
    matrixSize = size(matrix);
    level1BlockSize = matrixSize(1)/blockSize;
    level2BlockSize = matrixSize(2)/blockSize;
    totalDiagonalTerms = sum(diagonals);
    acTermsSize = level1BlockSize * level2BlockSize * totalDiagonalTerms; 
    acTerms = zeros(1, acTermsSize, 'double');
    keyCounter = 0;

    for i = 1:level1BlockSize
        for j = 1:level2BlockSize
            asc = true;
            for d = 1:length(diagonals)
                for t = 0:diagonals(d)-1
                    if asc
                        rowIndex = 1+t;
                        colIndex = diagonals(d)-t;
                    else
                        colIndex = 1+t;
                        rowIndex = diagonals(d)-t;
                    end

                    value = matrix((i - 1) * blockSize + rowIndex, (j - 1) * blockSize + colIndex);
                    key = (i - 1) * level2BlockSize  + j + keyCounter;
                    acTerms(key) = value;
                    keyCounter = keyCounter + 1;
                end
                asc = ~asc;
            end
        end
    end

    dFirst = floor(min(acTerms));
    dLast = ceil(max(acTerms) + power(10, -6));
    [dLevels, rLevels] = calcUniformIntervals(dFirst, dLast, 4);
    [acQuantized, acDequantized] = quanDequantArray(acTerms, dLevels, rLevels);
end

function [dequantizedImage] = dequantizeImage(matrix, blockSize, dcDequantized, acDequantized, diagonals)
    matrixSize = size(matrix);
    level1BlockSize = matrixSize(1)/blockSize;
    level2BlockSize = matrixSize(2)/blockSize;
    dequantizedImage = zeros(matrixSize(1), matrixSize(2));
    keyCounter = 0;

    for i = 1:level1BlockSize
        for j = 1:level2BlockSize
            dcKey = (i - 1) * level2BlockSize  + j;
            dequantizedImage((i - 1) * blockSize + 1, (j - 1) * blockSize + 1) = dcDequantized(dcKey);

            asc = true;
            for d = 1:length(diagonals)
                for t = 0:diagonals(d)-1
                    if asc
                        rowIndex = 1+t;
                        colIndex = diagonals(d)-t;
                    else
                        colIndex = 1+t;
                        rowIndex = diagonals(d)-t;
                    end
                    
                    acKey = (i - 1) * level2BlockSize  + j + keyCounter;
                    dequantizedImage((i - 1) * blockSize + rowIndex, (j - 1) * blockSize + colIndex) = acDequantized(acKey);
                    
                    keyCounter = keyCounter + 1;
                end
                asc = ~asc;
            end
        end
    end
end

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

function [quantized, dequantized] = quanDequantArray(array, dLevels, rLevels)
    arraySize = length(array);
    quantized = zeros(1, arraySize, 'double');
    dequantized = zeros(1, arraySize, 'double');

    for i = 1:arraySize
        dIndex = findLevelIndex(dLevels, array(i));
        quantized(i) = dIndex;
        dequantized(i) = rLevels(dIndex);
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
