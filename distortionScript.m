close all;
clear;
clc;

% adding folders and subfolders at path
addpath(genpath('bfilter2'));
addpath(genpath('Work Data'));
addpath(genpath('Help functions'));

work_imgs = dir('Work Data/*.jpg');
utils_dir = 'Utils';
utils_prefix = 'u_';

utilsPathPref = [utils_dir, '/', utils_prefix];
distortionMaskPath = [utilsPathPref, 'distortionMask.mat'];
distortionMaskImagePath = [utilsPathPref, 'distortionMaskImage.jpg'];
commonInfoPath = [utilsPathPref, 'commonInfo.mat'];

[~, ~] = mkdir(utils_dir);
photosCount = size(work_imgs, 1);
distortionOffset = [2350 1000];

% make distortion and save corresponding distortion mask
[~, ~, transparency] = imread('object.png');
distortion = imread('object.png');

I = imread([work_imgs(checkingImageIndex).folder, '/', work_imgs(checkingImageIndex).name]);

underlay = cutImageWithOffset(I, size(distortion), distortionOffset);
resultLay = insertWithTransparencyMask(underlay, distortion, transparency);
I = insertImageWithOffset(I, resultLay, distortionOffset);

distortionMask = zeros(size(I));
distortionMask = insertImageWithOffset(distortionMask, transparency, distortionOffset);
distortionMask(distortionMask > 0) = 1;
distortionMask = logical(sum(distortionMask, 3));

commonInfo.distortionOffset = distortionOffset;
commonInfo.photosCount = photosCount;
commonInfo.checkingImageIndex = checkingImageIndex;
commonInfo.filterType = filterType;
save(commonInfoPath, 'commonInfo');
save(distortionMaskPath, 'distortionMask');
imwrite(I, [utilsPathPref, work_imgs(checkingImageIndex).name]);
imwrite(distortionMask, distortionMaskImagePath);
% imshow(I);
