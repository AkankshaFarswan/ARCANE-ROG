clear
addpath(genpath('/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/sim_data_new/G_noisyData_nodoublets/GT_100x100_5_new/DT/')); 
%tic
Files=dir('/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/sim_data_new/G_noisyData_nodoublets/GT_100x100_5_new/DT/');
nhd=[];
impute_acc=[];
FPFN_ratio=[];
for kk=3:length(Files)
   FileNames=Files(kk).name
   mat_in=csvread(FileNames,0,0); 
   X=mat_in;
   ms=3; % ms represents missing data. In SNV data, if 3 represents missing, then ms=3; In CNV data, if -1 represents missing, then ms=-1.
   omega=find(X~=ms);
   omegaC=find(X==ms);
   [m,n]=size(X);
   if(m<=100)
       k=m/10;
   end
   if(m >100 && m<=500)
       k=m/25;
   end
   if(m>500 && m<=1000) 
       k=m/40;
   end
   alpha = 1/sqrt(max(m,n))*(1+3*length(omega)/(m*n));
   beta = 6/sqrt(max(m,n)); 
   mu=6/sqrt(max(m,n));
   %alpha=0.1521;
   %beta=0.1;
   %mu=0.5;
   %alphalist = [1/sqrt(max([m n]))];
   Y1 =zeros(m,n) ;
   Y2 =zeros(m,n) ;
   E = zeros(m,n);
   Z = X;
   distX = L2_distance_1(X,X);
   [distX1, idx] = sort(distX,2);
   count=10;
   [gamma] = cal_gamma(X,distX1,beta,k);
   count=count+1;
   disp([alpha,beta,mu])
   for i = 1:50
     D =  updateD(E,X,Y1,Y2,mu,Z,gamma);
     distX = L2_distance_1(D,D);
     [distX1, idx] = sort(distX,2);
     [gamma] = cal_gamma(D,distX1,beta,k);
     E = updateE_extend(D,E,X,Y1,mu,alpha,omega);
     S = updateS(X,distX1,idx,k,gamma,beta);
     S=(S+S')/2;
     L = diag(sum(S))-S;
     Z = updateZ(L,beta,mu,D,Y2);
     Y1 = Y1+mu*(D+E-X);
     Y2 = Y2+mu*(D-Z);
     mu=mu*1.1;
     if (norm(D+E-X,'inf')<1e-5) && (norm(D-Z,'inf')<1e-5)
         break
     end
   end
   L=(S+S')/2;
   AA1=int8(Z);
   %AA1=AA1';
   filedenoise = strcat(strtok(FileNames,'.'),'_denoised.mat')
   simi_mat = strcat(strtok(FileNames,'.'),'_S.mat')
   save(['/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/sim_data_new/G_noisyData_nodoublets/GT_100x100_5_new/DT_denoised_mRPCA_Smatrices/' simi_mat],'S')
   save(['/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/sim_data_new/G_noisyData_nodoublets/GT_100x100_5_new/DT_denoised_mRPCA/' filedenoise],'AA1')
   %csvwrite(['/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/matlab_and_R_scripts/scarletdata_denoisedMatrices/matfiles/'filedenoise],AA1)
   mat_denoised = double(AA1);
   temp1 = strtok(FileNames,'.')
   temp2 = strsplit(temp1,'_')
   strtemp = temp2{2}
   strfinal = strcat('GT_',strtemp,'.csv')
   mat_GD = csvread(fullfile('/mnt/storage/Akanksha/SingleCell/RobustClone/RobustClone-master/sim_data_new/G_noisyData_nodoublets/GT_100x100_5_new/GT/',strfinal),0,0);

   idx_true0 = find(mat_GD==0);
   idx_true1 = find(mat_GD==1);
   %%Error rate
   idx_unequal = find(mat_denoised ~= mat_GD);
   temp_nhd = size(idx_unequal,1)/(size(mat_GD,1)*size(mat_GD,2));
   nhd = [nhd;temp_nhd];
   
  
   %%FPs+FNs ratios of output GTM to input GTM
   in_true0 = mat_in(idx_true0);
   in_true1 = mat_in(idx_true1);
   FP_in = find(in_true0 == 1);
   FN_in = find(in_true1 == 0);
   out_true0 = mat_denoised(idx_true0);
   out_true1 = mat_denoised(idx_true1);
   FP_out = find(out_true0 == 1);
   FN_out = find(out_true1 == 0);
   temp_FPFNratio = (length(FP_out)+length(FN_out))/(length(FP_in)+length(FN_in));
   FPFN_ratio=[FPFN_ratio;temp_FPFNratio];
end