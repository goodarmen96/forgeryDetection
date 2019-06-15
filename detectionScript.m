close all;
clear;
clc;

% adding folders and subfolders at path
addpath(genpath('bfilter2'));
addpath(genpath('Work Data'));
addpath(genpath('Help functions'));

work_imgs = dir('Work Data/*.jpg');
work_noise_dir = 'Work Noise_';
fil_mah_dir = 'Mahalanobis Filtered';
utils_dir = 'Utils';
noised_prefix = 'noise_';
mah_prefix = 'mah_';
fil_mah_prefix = 'fil_mah_';
result_prefix = 'r_';
utils_prefix = 'u_';
workNoisePathPref = [work_noise_dir, '/', noised_prefix];
filMahPathPref = [fil_mah_dir, '/', fil_mah_prefix];
utilsPathPref = [utils_dir, '/', utils_prefix];
matPath = [utilsPathPref, '_mat.mat'];
dispersiaPath = [utilsPathPref, '_dispersia.mat'];
tresholdPath = [utilsPathPref, 'treshold.mat'];
resultPath = [utilsPathPref, 'result.mat'];
metricF1Path = [utilsPathPref, 'metricF1.mat'];
utilsNoisePathPref = [utilsPathPref, 'noise_'];
utilsMahPathPref = [utilsPathPref, 'mah_'];
distortionMaskPath = [utilsPathPref, 'distortionMask.mat'];
histImagePath = [utilsPathPref, 'histImage.jpg'];
resultImagePathPref = [utilsPathPref, 'resImage_'];
[~, ~] = mkdir(work_noise_dir);
[~, ~] = mkdir(fil_mah_dir);
photosCount = size(work_imgs, 1);

% photo sizes
photoSize = size(imread([work_imgs(1).name]));
photoHeight = photoSize(1);
photoWidth = photoSize(2);
photoSquare = photoHeight * photoWidth;

% size for mahalanobis distance filter
mahFilterSize = 21;
checkingImageIndex = photosCount;

% filterType: 1 - Bilateral, 2 - Gauss, 3 - Non-local means
filterType = 1;

% filter, getting noises
for i = 1:photosCount
    disp(['getting noises ', num2str(i), currentTime]);
    if i ~= checkingImageIndex 
        I = imread([work_imgs(i).folder, '/', work_imgs(i).name]);
    else
        I = imread([utilsPathPref, work_imgs(i).name]);
    end
    I = double(I)/255;
        
    % %{   
    % Apply filter to each image.
    flt_img = filtering(I, filterType);
    noise = I - flt_img;
    noise = (noise + 0.5)*255;
        
    if i ~= checkingImageIndex 
        save([workNoisePathPref, work_imgs(i).name, '.mat'], 'noise');
    else
        save([utilsNoisePathPref, work_imgs(i).name, '.mat'], 'noise');
    end
    
    % %}    
end

% compute mat of noises
mat = zeros(photoSize);
dispersia = mat;
for i = 1:photosCount
    if i == checkingImageIndex 
        continue;
    end
    disp(['mat ', num2str(i), currentTime]);
    lnoise = load([workNoisePathPref, work_imgs(i).name, '.mat']);
    noise = lnoise.noise;
    
    mat = mat + noise / (photosCount-1);
end

% compute dispersia of noises
for i = 1:photosCount    
    if i == checkingImageIndex 
        continue;
    end
    disp(['dispersia ', num2str(i), currentTime]);
    lnoise = load([workNoisePathPref, work_imgs(i).name, '.mat']);
    noise = lnoise.noise;
    
    dispersia = dispersia + (noise - mat).^2 / (photosCount-1);
end

subplot(1, 2, 1); imshow(uint8(linearContrast(mat))); title('noise mat')
subplot(1, 2, 2); imshow(uint8(dispersia)); title('noise dispersia')

save(matPath, 'mat');
save(dispersiaPath, 'dispersia');

% load each working noise and compute mahalanobis distance
loadMat = load(matPath, 'mat');
mat = loadMat.mat;
loadDisp = load(dispersiaPath, 'dispersia');
dispersia = loadDisp.dispersia;

mahFilter = ones([mahFilterSize mahFilterSize]);
centerIndex = uint8(mahFilterSize) / 2;
mahFilter(centerIndex, centerIndex) = sqrt(mahFilterSize);
mahFilter = mahFilter / sum(mahFilter(:)); 

for i = 1:photosCount
    disp(['mahalanobis ', num2str(i), currentTime]);
    if i ~= checkingImageIndex 
        loadWi = load([workNoisePathPref, work_imgs(i).name, '.mat']);
    else
        loadWi = load([utilsNoisePathPref, work_imgs(i).name, '.mat']);
    end
    Wi = loadWi.noise;
    mah = zeros(photoHeight, photoWidth);
    
    for n = 1:photoHeight
        for m = 1:photoWidth
            matV = [mat(n, m, 1) mat(n, m, 2) mat(n, m, 3)]';
            wiV = [Wi(n, m, 1) Wi(n, m, 2) Wi(n, m, 3)]';
            
            B = zeros(3, 3);
            B(1, 1) = dispersia(n, m, 1);
            B(2, 2) = dispersia(n, m, 2);
            B(3, 3) = dispersia(n, m, 3);

            dMah = (wiV - matV)'*B^-1*(wiV - matV);
            dMah = sqrt(dMah);

            mah(n, m) = dMah;
        end
    end

    if i ~= checkingImageIndex 
        mah  = imfilter(mah, mahFilter, 'circular');
        save([filMahPathPref, work_imgs(i).name, '.mat'], 'mah');
    else
        save([utilsMahPathPref, work_imgs(i).name, '.mat'], 'mah');
    end
end

% compute treshold using histogram from mahalanobis distances

% get histogram from all filtered mahalanobis distatnce
resultHistData = zeros([photoHeight, photoWidth * (photosCount-1)]);
histDataCount = 1;
for i = 1:photosCount
   if i == checkingImageIndex
       continue;
   end
   disp(['histogram ', num2str(i), currentTime]);
   loadIFilMah = load([filMahPathPref, work_imgs(i).name, '.mat']);
   iFilMah = loadIFilMah.mah;
   histRange = 1 + (histDataCount - 1) * photoSize(2) : histDataCount * photoSize(2);
   resultHistData(:, histRange) = iFilMah;
   
   histDataCount = histDataCount + 1;
end
figure; histogram(resultHistData); title('resultHistData');
saveas(gcf, histImagePath);
close figure 1;

% compute treshold
disp(['compute treshold', currentTime]);
treshold = tresholdFromDataPercent(resultHistData, 0.997);
save(tresholdPath, 'treshold');

% compute result binary of checked mahalanobis
loadTreshold = load(tresholdPath);
loadMah = load([utilsMahPathPref, work_imgs(checkingImageIndex).name, '.mat']);
treshold = loadTreshold.treshold;
iMah = loadMah.mah;  

result = zeros(size(iMah));
result(iMah > treshold) = 1;

CC = bwconncomp(result);
numOfPixels = cellfun(@numel, CC.PixelIdxList);

% remove areas which square less then 1/1000 of origin photo
numOfPixels(numOfPixels < (photoSquare / 1000)) = 0;
[~, indexes, ~] = find(numOfPixels > 0);
result = zeros(size(result));
for i = 1:length(indexes)
    result(CC.PixelIdxList{indexes(i)}) = 1;
end

% convex hull
result = bwconvhull(result, 'objects');
result = logical(result);
% figure; imshow(result);
imwrite(result, [resultImagePathPref, work_imgs(checkingImageIndex).name]);
save(resultPath, 'result');

% compute F1 metrics
loadR = load(resultPath);
loadMask = load(distortionMaskPath);
result = loadR.result;
distortionMask = loadMask.distortionMask;

TP = and(result, distortionMask);
FN = and(~result, distortionMask);
FP = and(result, ~distortionMask);

TP = sum(TP(:));
FN = sum(FN(:));
FP = sum(FP(:));

Precision = TP / (TP + FP);
Recall = TP / (TP + FN);

metricF1 = 2 * Precision * Recall / (Precision + Recall);
save(metricF1Path, 'metricF1');

disp(['end', currentTime]);

% some helps functions

% input - double [0 1], output - double [0 1]
function [filtered] = filtering(I, type)
if type == 1
    % Set bilateral filter parameters.
    w     = 5;       % bilateral filter half-width
    sigma = [3 0.1]; % bilateral filter standard deviations

    filtered = bfilter2(I, w, sigma);
elseif type == 2
    filtered = imgaussfilt(I); 
elseif type == 3
    filtered = imnlmfilt(I);
end

end
