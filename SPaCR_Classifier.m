function class = SPaCR_Classifier(DataTrains, CTrain, DataTests, lambda, bata, seg)

DataTrain = DataTrains(:,3:end);
DataTest  = DataTests(:,3:end); 

DDT = DataTrain*DataTrain';

numClass = length(CTrain);
[m Nt]= size(DataTest);
% information transition from pixel to superpixel, including spectral
% reflection and coordinate
superlabel = seg.labels;
S_C = seg.X_c';
Datasample = [DataTrains(:,1:2);DataTests(:,1:2)];
LinearInd = superlabel(sub2ind(size(superlabel),Datasample(:, 1),Datasample(:, 2)));
SuperInfo = S_C(LinearInd,:);

SuperDataTrain = SuperInfo(1:size(DataTrains,1),1:end-3);
SuperDataTest = SuperInfo(size(DataTrains,1)+1:end,1:end-3);

SuperTrainInd = SuperInfo(1:size(DataTrains,1),end-2:end-1);
SuperTestInd = SuperInfo(size(DataTrains,1)+1:end,end-2:end-1);


for j = 1: m
    if mod(j,round(m/20))==0
        fprintf('# -SPaCR is interpreting hyperspectral data...\n');        
    end
    xy = SuperTestInd(j,:); 
    XY = SuperTrainInd; 
    Spa_dist = sum((abs(XY' - repmat(xy', [1 size(XY,1)]))).^2);
    Spa_dist = Spa_dist./max(Spa_dist);
    xf = SuperDataTest(j,:);
    XF = SuperDataTrain;
    Spe_dist = pdist_le(XF,xf,'SAM');
    Spe_dist = Spe_dist./max(Spe_dist);
    Spe_dist = Spe_dist';
    gamma = 10^-6;
    spa_spe_dist = sqrt(Spa_dist.^2 + gamma*Spe_dist.^2);
    spa_spe_dist = spa_spe_dist./max(spa_spe_dist);
    D = diag(bata.*spa_spe_dist);
    
    Y = DataTest(j, :); % 1 x dim
    norms = sum((DataTrain' - repmat(Y', [1 size(DataTrain,1)])).^2);
    % norms = ones(size(DataTrain,1), 1);
    G = diag(lambda.*norms);
    weights = (DDT +  G + D)\(DataTrain*Y');
    
    a = 0;
    for i = 1: numClass 
        % Obtain Multihypothesis from training data
        HX = DataTrain((a+1): (CTrain(i)+a), :); % sam x dim
        HW = weights((a+1): (CTrain(i)+a));
        a = CTrain(i) + a;
        Y_hat = HW'*HX;
        
        Dist_Y(j, i) = norm(Y - Y_hat); 
    end
   Dist_Y(j, :) = Dist_Y(j, :)./sum(Dist_Y(j, :));
end
[~, class] = min(Dist_Y'); 
