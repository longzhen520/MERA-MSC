clear all; 
close all; 
clc;

addpath('ClusteringMeasure', 'mylib', 'datasets');
prompt = " Please choose the dataset name:'Yale' or 'MSRC' or 'EYaleB' or 'NH' or 'BDGP' \n";
dataname = input(prompt);
data=dataname;
switch data
    case 'Yale'      
         load('yale.mat');
          rX= [11 15 15 11 3];
          R=6;
          lambda=0.001;
          
    case 'MSRC'      
         load('MSRC.mat');
          rX= [15 14 15 14 5];
          R=2;
          lambda=0.001;
          
     case 'EYaleB'      
         load('EYaleB10_mtv.mat');
          rX= [32 20 20 32 3];
          R=10;
          lambda=1;
          
     case 'NH'  
         %% needs to download Notting Hill datasets
         load('NH.mat');
          rX= [233 10 233 10 3];
          R=15;
          lambda=0.0001;
          X{1}=X1(:,2:2:end);
          X{2}=X2(:,2:2:end);
          X{3}=X3(:,2:2:end);
          gt=gt(2:2:end);
      case 'BDGP'      
          load('BDGP_4view.mat')
          rX= [50 50 50 50 4];
          R=10;
          lambda=0.0002;
          gt=labels;
          V=length(X);
        for v=1:V
           X{v}=X{v}';
        end
    
          
          
end



    
%% Note: each column is an sample (same as in LRR)
V=length(X);
cls_num = length(unique(gt));
%data preparation...

for v=1:V
%     X{v}=fea{v}';
    [X{v}]=NormalizeData(X{v});
     %X{v} = zscore(X{v},1);
end
% Initialize...
%% need to tune

%% parameter setting
K = length(X); N = size(X{1},2); %sample number
paras_mera.R{1}=[R,R];
%    paras_mera.R{1}=[10,10];
paras_mera.lambda=lambda;
%paras_mera.M=floor(A(a)*N);
paras_mera.rX   = rX;
        
% -------------------0。5------------------- clustering 
   tic;
   [S,mera]=MERA_MSC(X,paras_mera);
   Time=toc;
    for q=1:10
       C=SpectralClustering(S, cls_num);
  % C = kmeans(U,numClust,'EmptyAction','drop');
        [Fi(q),Pi(q),Ri(q)] = compute_f(gt,C);
        [A1 nmi1(q) avgenti(q)] = compute_nmi(gt,C);    
        ACCi(q) = Accuracy(C,double(gt));
        if (min(gt)==0)
            [ARi(q),RIi(q),MIi(q),HIi(q)]=RandIndex(gt+1,C);
        else
            [ARi(q),RIi(q),MIi(q),HIi(q)]=RandIndex(gt,C);
        end       
     end
     
        F= mean(Fi); VF= std(Fi);
        P= mean(Pi); VP= std(Pi);
        Rx= mean(Ri); VR= std(Ri);
        nmi= mean(nmi1); Vnmi= std(nmi1);
        avgent= mean(avgenti); avgent= std(avgenti);
        AR= mean(ARi); VAR= std(ARi);
        ACC=mean(ACCi); VACC=std(ACCi);  
        
        
        
%          fprintf('S\n'); 
     fprintf('F-score: %.3f(%.3f)\n', F,VF);
    fprintf('P: %.3f(%.3f)\n', P,VP);    
    fprintf('R: %.3f(%.3f)\n', Rx,VR);
    fprintf('nmi:%.3f(%.3f)\n', nmi,Vnmi);
    fprintf('AR: %.3f(%.3f)\n', AR,VAR);
    fprintf('ACC: %.3f(%.3f)\n', ACC,VACC);  

