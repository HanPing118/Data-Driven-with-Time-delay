clear;clc;tic;
a=load('E:\BaiduSyncdisk\aaa临时\PAPER1\CHENGXU\Paper_1\Data_Driven_delay\Sample_delay2_2.txt');
x=a(:,1);y1=a(:,2);y2=a(:,3);wbin=a(:,4);  %漂移和扩散系数和权重

G=100;    %  对应fortran的nbin
L = length(find(isnan(y2)));     % 数据预处理
wbin(isnan(y2)) = [];
x(isnan(y2)) = [];
y1(isnan(y2)) = [];
y2(isnan(y2)) = [];

% figure(10);
% tau=0.1;DD=0.05;xx=linspace(0,3,100);yy1=((1-2*tau*DD).*xx-xx.^2)./(1+tau.*xx);
% plot(x,y1,'-','LineWidth',3,'MarkerSize',6);hold on
% plot(xx,yy1,'-','LineWidth',3,'MarkerSize',6);hold on
% % figure(11);
% xx=linspace(0,3,100);yy2=xx.^2;
% plot(x,y2,'-','LineWidth',3,'MarkerSize',6);hold on
% plot(xx,yy2,'-','LineWidth',3,'MarkerSize',6);hold on

wbin=diag(wbin);
M=4;Theta=zeros(G-L,M);
for i=1:M
    Theta(:,i)=x.^(i-1);
end
% A=ones(G-L,1);M=6;
% Theta=[A x x.^2 x.^3 x.^4 x.^5];

Theta1=wbin*Theta;
Theta2=wbin*Theta;
y1=wbin*y1;
y2=wbin*y2;
Xi1=pinv((Theta1'*Theta1))*Theta1'*y1;             % 标准最小二乘法
Xi2=pinv((Theta2'*Theta2))*Theta2'*y2;
%%*** CV+SSR 交叉验证 Cross validation and Stepwise Sparse regressor
K=10;
E_Vec1=zeros(M,K);E_Vec2=zeros(M,K);
Delta_Vec1=zeros(M,1);Delta_Vec2=zeros(M,1);
Ave_train_Xi1=zeros(K*M,M);Ave_train_Xi2=zeros(K*M,M);
Ave_Xi1=zeros(M,M);Ave_Xi2=zeros(M,M);
%将数据样本随机分割为K=10部分,K折交叉
indices1=crossvalind('Kfold',Theta1(1:G-L,M),K); 
indices2=crossvalind('Kfold',Theta2(1:G-L,M),K);
fid1=fopen('Position_Drift.txt','w');
for q = 1:M                                    %q是参数为0的个数
    G1=0;
    fprintf(fid1, '%25s  %d\r\n', 'the number of coefficient',M-q);
    for k = 1:K
        test = (indices1 == k);    % 获取第i份测试数据的索引逻辑值
        train = ~test;   % 取反，获取第i份训练数据的索引逻辑值  

        test_data_Theta1 = Theta1(test,:);   %1份测试，9份训练
        test_label_y1 = y1(test);
        
        train_data_Theta1 = Theta1(train,:);
        train_label_y1 = y1(train);
        % SSR
        train_Xi1= pinv((train_data_Theta1'*train_data_Theta1))*train_data_Theta1'*train_label_y1;
        for j=1:q
            smallinds1=find(abs(train_Xi1)==min(abs(train_Xi1)));
            smallinds1=smallinds1(1);
            
            train_data_Theta1(:,smallinds1)=[];   %making the smallest one be 0
            test_data_Theta1(:,smallinds1)=[];
            % Regress dynamics onto remaining terms to find sparse Xi
            train_Xi1= pinv((train_data_Theta1'*train_data_Theta1))*train_data_Theta1'*train_label_y1;  
            if j == q
                fprintf(fid1, '%d\r\n',smallinds1);
            else
                fprintf(fid1, '%d  ',smallinds1);
            end
        end
        %记录中间参数
        [m1,n1]=size(train_Xi1);
        Ave_train_Xi1((q-1)*K+k,1:m1)=train_Xi1;
       % 交叉验证得分
        E_Vec1(q,k)=norm(test_label_y1-test_data_Theta1*train_Xi1,2);
        G1=G1+E_Vec1(q,k); 
    end
    Ave_Xi1(q,1:m1)=mean(Ave_train_Xi1((q-1)*K+1:(q-1)*K+K,1:m1));
    Delta_Vec1(q)=G1/K;
end 
%将数据样本随机分割为K=10部分,K折交叉
fid2=fopen('Position_Diffusion.txt','w');
for q = 1:M                                    %q是参数为0的个数
    G2=0;
    fprintf(fid2, '%25s  %d\r\n', 'the number of coefficient',M-q);
    for k = 1:K
        test = (indices2 == k);    % 获取第i份测试数据的索引逻辑值
        train = ~test;   % 取反，获取第i份训练数据的索引逻辑值  

        test_data_Theta2 = Theta2(test,:);
        test_label_y2 = y2(test);
        
        train_data_Theta2 = Theta2(train,:);
        train_label_y2 = y2(train);
        % SSR
        train_Xi2= pinv((train_data_Theta2'*train_data_Theta2))*train_data_Theta2'*train_label_y2;
        for j=1:q
            smallinds2=find(abs(train_Xi2)==min(abs(train_Xi2)));
            smallinds2=smallinds2(1);
            
            train_data_Theta2(:,smallinds2)=[];
            test_data_Theta2(:,smallinds2)=[];
            % Regress dynamics onto remaining terms to find sparse Xi
            train_Xi2= pinv((train_data_Theta2'*train_data_Theta2))*train_data_Theta2'*train_label_y2;
            if j == q
                fprintf(fid2, '%d\r\n',smallinds2);
            else
                fprintf(fid2, '%d  ',smallinds2);
            end
        end
        %记录中间参数
        [m2,n2]=size(train_Xi2);
        Ave_train_Xi2((q-1)*K+k,1:m2)=train_Xi2;
        % 交叉验证得分
        E_Vec2(q,k)=norm(test_label_y2-test_data_Theta2*train_Xi2,2);
        G2=G2+E_Vec2(q,k); 
    end
    Ave_Xi2(q,1:m2)=mean(Ave_train_Xi2((q-1)*K+1:(q-1)*K+K,1:m2));
    Delta_Vec2(q)=G2/K;
end 
save apr_smooth_test1.mat Delta_Vec1 Delta_Vec2 M Ave_Xi1 Ave_Xi2 
% matlab中的flipud函数实现矩阵的上下翻转。flipud(X)实现了矩阵X的上下翻转。 
figure;
plot(1:M-1, flipud(Delta_Vec1(1:M-1)), 'b-*','LineWidth',2,'MarkerSize',6);hold on
% figure;
plot(1:M-1, flipud(Delta_Vec2(1:M-1)), 'r-*','LineWidth',2,'MarkerSize',6);hold on

% set(gca,'YTick');
% set(gca,'YTicklabel',{'10^-^4','10^-^3','10^-^2','10^-^1','10^0','10^1','10^2'},'Fontname', 'Times New Roman','FontSize',12)
legend('Drift term','Diffusion term')
xlabel('\itn')
ylabel('\it\delta(\itn)')
title 'Cross validation score'

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M=8;K=10;
% Coef_1=[1  6  4  2  4
% 7  1  4  2  4
% 7  1  4  2  4
% 5  1  5  2  4
% 1  6  4  2  4
% 7  3  4  1  4
% 7  1  4  2  4
% 7  1  4  2  4
% 7  1  4  2  4
% 7  3  4  1  4];
% Coef_2=[2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2
% 2  7  5  3  3  3  2];
% for i=1:K
%     F1=Find_basis(Coef_1(i,:),M);
%     F2=Find_basis(Coef_2(i,:),M);
%     Final1{i}=F1;
%     Final2{i}=F2;
% end
% celldisp(Final1)
% celldisp(Final2)
% save basis1.mat  Final1
% save basis2.mat  Final2
%     
% function F=Find_basis(Coef,M)
%     basis={'1',' x','x.^2','x.^3','x.^4','x.^5','x.^6','x.^7'};
%     pos=1:1:M;
%     for i=1:1:length(Coef)
%         pos(Coef(i))=[];
%     end
%     for j=1:1:length(pos)
%         F{j}=basis{pos(j)};
%     end    
% end 