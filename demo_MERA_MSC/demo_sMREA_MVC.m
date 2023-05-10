clear all; 
close all; 
clc;

addpath('ClusteringMeasure', 'mylib', 'datasets');
prompt = " Please choose the dataset name:'ALOI' or 'Caltech' or 'CCV' \n";
dataname = input(prompt);
data=dataname;
switch data
    case 'ALOI'      
         load('ALOI_100.mat');
          rX=[108 100 108 100 4];
           V=length(fea);
               R=4;
          lambda=0.02;
               M=10; %% anchor numbers
        for v=1:V
           X{v}=fea{v}';
        end
        
    
          
    case 'Caltech'      
         load('Caltech101-all_fea.mat');         
         rX=[114 80 114 80 5];
          R=8;
          lambda=0.1;
          M=9;
      V=length(X);
        for v=1:V
           X{v}=X{v}(1:end-24,:)';
        end
        gt=Y(1:end-24,1);
        
     case 'CCV'      
         load('CCV.mat');
          rX= [67 100 67 100 3];
          R=2;
          lambda=0.5;
          M=6;
      V=length(X);
        for v=1:V
           X{v}=X{v}(1:end-73,:)';
        end
        gt=Y(1:end-73,1);
          
     
    
          
          
end

%% load data

    
%% Note: each column is an sample (same as in LRR)
V=length(X);
cls_num = length(unique(gt));
%data preparation...

for v=1:V
    [X{v}]=NormalizeData(X{v});
end
% Initialize..


   
   
      
    

%% parameter setting
K = length(X); N = size(X{1},2); %sample number
paras_mera.R{1}=[R,R];
paras_mera.lambda=lambda;
paras_mera.M   = M;
paras_mera.rX   = [M,rX(3),rX(4),rX(5)];

        
% -------------------0。5------------------- clustering 
   tic;
   [S,mera]=scale_MERA_MVC(X,paras_mera);
   Time=toc;
   
    [U,~,~]=svd(S','econ');
   
    for q=1:10
       C=kmeans(U, cls_num);
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
   fprintf('CPU Time: %.3f\n', Time);  

