clc;
clear;
addpath(genpath('./utils/'));
addpath(genpath('./codes/'));

result_URL = './results/';
% if ~isfolder(result_URL)
%     mkdir(result_URL);
% end
dataset_name = {'nus-wide'};
%   rng('default');
%% load dataset
 for db = 1:length(dataset_name)
  dataset = dataset_name{db};
  load(['./datasets/',dataset,'.mat']);  
%% parameter setting

     I_tr=XDatabase;
     I_te=XTest;
     T_tr=YDatabase;
     T_te=YTest;
     L_tr=databaseL;
     L_te=testL;
%% 3  
 nbits = [4];
%  mbits = ceil(nbits*0.25);
%  mbits = [0:8];
 beta = [0.1];
 lambda = [10];
 muta = [10];
 theta = [0.001];
 maxItr =  5;
 turn = 3; %turns
 func = 'linear';
%    [d,c]=size(L_tr); 
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
                    CBHparam.nAnchors = 2000; 
                   [XKTrain,XKTest] = Kernelize(I_tr,I_te,CBHparam.nAnchors); 
                   [YKTrain,YKTest] = Kernelize(T_tr, T_te,CBHparam.nAnchors);
                    XTrain = XKTrain';
                    YTrain = YKTrain';
                    LTrain = L_tr;
                    LTest = L_te;
                    CBHparam.lambda = lambda(j);
                    CBHparam.muta = muta(k);
                    CBHparam.theta= theta(h);
                    CBHparam.nbits = nbits(bi);
                    CBHparam.mbits = ceil(nbits(bi)*0.20);
                    CBHparam.beta = beta(r);
                    CBHparam.maxItr = maxItr(b);
                    CBHparam.func = func;
                    [B] = CBH(XTrain,YTrain,LTrain',CBHparam);
                    eva_info_ = evaluate_CBH(XTrain',YTrain',LTrain,XKTest,YKTest,LTest,CBHparam,B);
                    eva_info_.Image_VS_Text_MAP;
                    eva_info_.Text_VS_Image_MAP;
                    result.bits = nbits(bi);
                    result.muta = muta(k);
                    arry(l,1) = nbits(bi);
                    arry(l,3) = lambda(j);
                    arry(l,4) = muta(k);
                    arry(l,6) = eva_info_.Image_VS_Text_MAP;
                    arry(l,7) = eva_info_.Text_VS_Image_MAP;
                    arry(l,8) = eva_info_.Image_VS_Text_MAP + eva_info_.Text_VS_Image_MAP; 
                    l=l+1;
                     map(v,1)=eva_info_.Image_VS_Text_MAP;
                     map(v,2)=eva_info_.Text_VS_Image_MAP;
                     top(v,1) = eva_info_.I2Ttop;
                     top(v,2) = eva_info_.T2Itop;
                  end
             MAP.i2t(b) = mean(map( : , 1));
             MAP.t2i(b) = mean(map( : , 2));
             MAP.Mean(b) = (MAP.i2t(b)+MAP.t2i(b))/2;
             TOP.i2t(b) = mean(top( : , 1));
             TOP.t2i(b) = mean(top( : , 2));
             TOP.Mean(b) = (TOP.i2t(b)+TOP.t2i(b))/2;
             Train_time = mean(t2);
             Time = mean(t1);
             fprintf('%d bits average map over %d runs for ImageQueryForText: %.4f\n, top@100 is %.4f\n',nbits(bi), turn, mean(map( : , 1)),mean(top( : , 1)) );
             fprintf('%d bits average map over %d runs for TextQueryForImage:  %.4f\n, top@100 is %.4f\n',nbits(bi), turn, mean(map( : , 2)),mean(top( : , 2)));
               end 
             end
           end 
        end       
    end


end

 end


