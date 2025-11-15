clc;
clear;
addpath(genpath('./utils/'));
addpath(genpath('./codes/'));
result_URL = './results/';
dataset_name = {'mirflickr25k'};

%% load dataset
 for db = 1:length(dataset_name)
  dataset = dataset_name{db};
  load(['./datasets/',dataset,'.mat']);  
%% parameter setting 
 nbits = [4];
 lambda = [10];
 muta = [10];
 beta = [0.1];
 theta = [0.001];
 maxItr = 5;
 turn = 10; %turns 
 func = 'linear';
   [d,c]=size(L_tr); 
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
                    XTest = XKTest;
                    YTest = YKTest;
                    LTest = L_te;
                    
                    CBHparam.lambda = lambda(j);
                    CBHparam.muta = muta(k);
                    CBHparam.theta= theta(h);
                    CBHparam.nbits = nbits(bi);
%                   CBHparam.mbits = mbits(b);
                    CBHparam.mbits = ceil(CBHparam.nbits*0.2);
                    CBHparam.beta = beta(r);
                    CBHparam.maxItr = maxItr(b);
                    CBHparam.func = func;
                    [B] = CBH(XTrain,YTrain,LTrain',CBHparam);
                    eva_info_ = evaluate_CBHH(XTrain',YTrain',LTrain,XTest,YTest,LTest,CBHparam,B);
                    eva_info_.Image_VS_Text_MAP;
                    eva_info_.Text_VS_Image_MAP;
                    result.bits = nbits(bi);
                    result.muta = muta(k);
                    map(v,1)=eva_info_.Image_VS_Text_MAP;
                    map(v,2)=eva_info_.Text_VS_Image_MAP;
                    top(v,1) = eva_info_.I2Ttop;
                    top(v,2) = eva_info_.T2Itop;
                  end
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


