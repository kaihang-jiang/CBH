function final_B = CBH(XKTrain,YKTrain,LTrain,param)
  %%  parameters
     
     lambda = param.lambda;
     muta = param.muta; 
     nbits = param.nbits;
     mbits = param.mbits;
     beta = param.beta*nbits/4;
     cbits = (nbits-mbits*2);
     sbits = cbits+mbits;
     f = 13/24;
%            rng('default');
      [c,n] = size(LTrain);
      [d1,~] = size(XKTrain);
      [d2,~] = size(YKTrain);
  %%  initial
       B = sgn(0.6*randn(nbits,n));
       C = B(1:cbits,:);
       C1 = B(cbits+1:sbits,:);
       C2 = B(sbits+1:nbits,:);
      LTrain_ = LTrain;
      LTrain_(LTrain_==0)=-1;
      LTrain_(LTrain_==1)=0;
      LTrain_(LTrain_==-1)=1;
      F = LTrain;
      for i = 1:n
      if norm(LTrain(:,i),1) == 0
      F(:,i) = 0;
      F_(:,i) = 1;
      else
      F(:,i) = LTrain(:,i)/norm(LTrain(:,i),1);
      F_(:,i) = LTrain_(:,i)/norm(LTrain_(:,i),1);
      end
      end  
      G = NormalizeFea(LTrain,0);
      G1  = [F;LTrain;LTrain;F_];
      G2 = [LTrain;F;-F_;-LTrain]/2; 
      XKTrain=NormalizeFea(XKTrain,0);
      YKTrain=NormalizeFea(YKTrain,0);
    for i = 1:param.maxItr 
     %% update Wt 
       temp_1 = ([C;C1]*XKTrain');
       [U_1,~,V_1] = svd(temp_1,'econ');
       W_1 = U_1*V_1';
       temp_2 = ([C;C2]*YKTrain');
       [U_2,~,V_2] = svd(temp_2,'econ');
       W_2 = U_2*V_2';
%        
     %% update V 
        V1 = (sbits*([C;C1]*G1')*G2)';
        Temp = V1'*V1-1/(n)*(V1'*ones(n,1)*(ones(1,n)*V1));
        [~,Lmd,RR] = svd(Temp);
        idx = (diag(Lmd)>1e-5);
        R = RR(:,idx); R_ = orth(RR(:,~idx));
        P = (V1-1/(n)*ones(n,1)*(ones(1,n)*V1)) *  (R / (sqrt(Lmd(idx,idx))));
        P_ = orth(randn(n,sbits-length(find(idx==1))));
        U1 = sqrt(beta*n)*[P P_]*[R R_]';
        V1 = (U1)'; 
        clear idx RR Lmd Temp
        
        V2 = (sbits*([C;C2]*G1')*G2)';
        Temp = V2'*V2-1/(n)*(V2'*ones(n,1)*(ones(1,n)*V2));
        [~,Lmd,RR] = svd(Temp);
        idx = (diag(Lmd)>1e-5);
        R = RR(:,idx); R_ = orth(RR(:,~idx));
        P = (V2-1/(n)*ones(n,1)*(ones(1,n)*V2)) *  (R / (sqrt(Lmd(idx,idx))));
        P_ = orth(randn(n,sbits-length(find(idx==1))));
        U2 = sqrt(beta*n)*[P P_]*[R R_]';
        V2 = (U2)';  
        clear idx RR Lmd Temp
        
     %% update C
        wW1 = W_1(1:cbits,:);
        bW1 = W_1(cbits+1:end,:);
        wW2 = W_2(1:cbits,:);
        bW2 = W_2(cbits+1:end,:);
        wV1 = V1(1:cbits,:);
        bV1 = V1(cbits+1:end,:);
        wV2 = V2(1:cbits,:);
        bV2 = V2(cbits+1:end,:);
        wB = B(1:cbits,:);
        C = sgn(((muta*beta*2*n+lambda*2+1)*eye(cbits))\...
            (lambda*(wW1*XKTrain+wW2*YKTrain)+muta*((sbits*f*(wV1*G1')*G2-wV1*bV1'*C1)+(sbits*f*(wV2*G1')*G2-wV2*bV2'*C2))+wB));
        
      %% update Ct
        B1 = B(cbits+1:sbits,:);
        B2 = B(sbits+1:nbits,:);
        C1 = sgn(((muta*beta*n+lambda*2+1)*eye(mbits))\...
            (lambda*(bW1*XKTrain-C2)+ muta*(sbits*f*(bV1*G1')*G2-bV1*wV1'*C)+B1));
        C2 = sgn(((muta*beta*n+lambda*2+1)*eye(mbits))\...
            (lambda*(bW2*YKTrain-C1)+ muta*(sbits*f*(bV2*G1')*G2-bV2*wV2'*C)+B2));

        %% update B    
        
        B = ([C; C1; C2]);     
    end  
     final_B = B;

end
