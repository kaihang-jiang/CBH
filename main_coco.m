clc;
clear;
addpath(genpath('./utils/'));
addpath(genpath('./codes/'));

result_URL = './results/';
% if ~isfolder(result_URL)
%     mkdir(result_URL);
% end
dataset_name = {'coco'};
%   rng('default');
%% load dataset
 for db = 1:length(dataset_name)
  dataset = dataset_name{db};
  load(['./datasets/',dataset,'.mat']);  
%% parameter setting 
 nbits = [4];
 lambda = [10];
 muta = [10];
 beta = [0.05];
 theta = [0.001];
 maxItr = 5;
 turn = 5; %turns 
 func = 'linear';
 [d,c]=size(L_tr); 
l = 1; %excel writing parameter
 %% start
for bi = 1:length(nbits)
    for j=1:length(lambda)
        for b = 1:length(maxItr)
           for h = 1:length(theta)
             for k = 1:length(muta)
               for r=1:length(beta)
                  for v = 1:turn
                    fprintf('Data preparing...\n\n'); 
                    XTrain = I_tr';
                    YTrain = T_tr';
                    LTrain = L_tr;
                    XTest = I_te;
                    YTest = T_te;
                    LTest = L_te;
                    
                    BHSHparam.lambda = lambda(j);
                    BHSHparam.muta = muta(k);
                    BHSHparam.theta= theta(h);
                    BHSHparam.nbits = nbits(bi);
                    BHSHparam.mbits = ceil(BHSHparam.nbits*0.2);
                    BHSHparam.beta = beta(r);
                    BHSHparam.maxItr = maxItr(b);
                    BHSHparam.func = func;
                    [B] = CBH(XTrain,YTrain,LTrain',BHSHparam);
                    eva_info_ = evaluate_CBH(XTrain',YTrain',LTrain,XTest,YTest,LTest,BHSHparam,B);
                    Bx = eva_info_.BX;  By = eva_info_.BY;
                    eva_info_.Image_VS_Text_MAP;
                    eva_info_.Text_VS_Image_MAP;
                    result.bits = nbits(bi);
                    result.muta = muta(k);
%                     result.theta =theta(r);
                    arry(l,1) = nbits(bi);
                    arry(l,3) = lambda(j);
                    arry(l,4) = muta(k);
%                     arry(l,5) = theta(r);
                    arry(l,6) = eva_info_.Image_VS_Text_MAP;
                    arry(l,7) = eva_info_.Text_VS_Image_MAP;
                    arry(l,8) = eva_info_.Image_VS_Text_MAP + eva_info_.Text_VS_Image_MAP; 
                    l=l+1;
                     map(v,1)=eva_info_.Image_VS_Text_MAP;
                     map(v,2)=eva_info_.Text_VS_Image_MAP;
                     top(v,1) = eva_info_.I2Ttop;
                     top(v,2) = eva_info_.T2Itop;
                  end
              mmapi2t(j,k) = mean(map( : , 1));
              mmapt2i(j,k) = mean(map( : , 2));
             fprintf('%d bits average map over %d runs for ImageQueryForText: %.4f\n, top@100 is %.4f\n',nbits(bi), turn, mean(map( : , 1)),mean(top( : , 1)) );
             fprintf('%d bits average map over %d runs for TextQueryForImage:  %.4f\n, top@100 is %.4f\n',nbits(bi), turn, mean(map( : , 2)),mean(top( : , 2)));
               end 
             end
           end
             MAP.i2t(b) = mean(map( : , 1));
             MAP.t2i(b) = mean(map( : , 2));
             MAP.Mean(b) = (MAP.i2t(b)+MAP.t2i(b))/2;
             TOP.i2t(b) = mean(top( : , 1));
             TOP.t2i(b) = mean(top( : , 2));
             TOP.Mean(b) = (TOP.i2t(b)+TOP.t2i(b))/2;
        end       
    end
end

 end

