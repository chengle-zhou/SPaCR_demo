% =========================================================================
% Spatial Peak-Aware Collaborative Representation for Hyperspectral Imagery Classification
% MATLAB demo:
%
% Reference
% [1] C. Zhou, B. Tu, Q. Ren, and S. Chen, ¡°Spatial Peak-Aware CollaborativeRepresentation
%     for Hyperspectral Remote Sensing Image Classification,¡± IEEE Geoscience and Remote Sensing Letters,
%     2021, to be published,doi:10.1109/LGRS.2021.3083416
%
% For any questions, email me by chengle_zhou@foxmail.com
%=========================================================================

clear all; close all; clc
addpath('utilities')
addpath('.\data');
load Indian_pines_corrected;load Indian_pines_gt;

data = img;
gth = indian_pines_gt;
no_class = max(gth(:));
[row, col, band] = size(img);
row_data_2d = reshape(data,row*col,band);
row_data_2d = normcol(row_data_2d')'; % L2 norm normalization; usage: m(dim) * n(num)

% 1. SLIC image segmentation
% Sw - the averange width of superpixels; {3-11} 
% P  - the number of superpixels;
% Ws - Trade-off coefficient between spatial and spectral distances;{0.5}
if(~exist('superpixelinfo.mat'))
    Sw = 3; P = round(row*col/Sw^2); Ws = 0.2;
    seg = slic_HSI_le(data, P, Ws);
    save superpixelinfo seg
else
    load superpixelinfo
end
% 2. randomly divide the dataset to training and test samples
CTrain = [3,13,9,3,5,6,3,5,3,8,25,6,3,11,4,3]; % the rate of training is set to 1.0% of ground thruth
indexes = train_random_select(GroundT(2,:),CTrain); 
train_SL = GroundT(:,indexes);
test_SL = GroundT;
test_SL(:,indexes) = [];
[Cx Cy] = ind2sub([row col],1:row*col);
Coor_xy = [Cx',Cy'];
DataTrain = [Coor_xy(train_SL(1,:),:),row_data_2d(train_SL(1,:),:)];
DataTest = [Coor_xy(test_SL(1,:),:),row_data_2d(test_SL(1,:),:)];

% 3. parameter settings
lambda = 0.2;
bata  = 0.4;

% 4. SPaCR Classification
tic
test_result = SPaCR_Classifier(DataTrain, CTrain, DataTest, lambda, bata, seg);
SdCR_time = toc;
[OA,kappa,AA,CA] = calcError(test_SL(2,:)-1, test_result-1, 1:no_class)

% 5. Display the result of SPaCR
SPaCRresults = gth;
SPaCRresults(test_SL(1,:)) = test_result;
SPaCR_Map=label2color(SPaCRresults,'india');
figure,imshow(SPaCR_Map);

