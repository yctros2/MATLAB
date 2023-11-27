function[A,N,n,a,b,m,flag,mA,nA]=simplex(A,m,b,c)%A为初矩阵，N为A的最后一行，n为初始进基变量，a为出基，b为进基，m为最后基变量，flag为cj标志，mA，nA为A的行数，列数%此处求的最小极值；
% 输入:
%   A - 初始系数矩阵
%   m - 基变量
%   b - 约束条件右侧向量
%   c - 目标函数系数向量
% 输出:
%   A, N, n, a, b, m, flag, mA, nA - 单纯形表相关参数
A=[3 -1 2;-2 4 0;-4 3 8];
b=[7;12;12];
c=[1 -3 -1];
A= sym(A);
b = sym(b);
c = sym(c);

 %化成标准型；
[a1 a2]=size(A);    %取A的行数a1,列数a2；
B=eye(a1);        %取与行数相等的单位阵；
A1=[A B];         %标准单纯形自变量系数矩阵；
m=[(a2+1):(a1+a2)]; %列出基变量的下标号；
mm=[1:a2];         %列出非基变量下标；
m1=zeros(1,((a1+a2)-size(c,2)));%给c补充0向量；
c1=[c m1];                %生成完整的c向量；
AA1=[A1 b];
flag=1;
k=0;
rl=0;
r=1;
while flag==1
b=AA1(:,a1+a2+1);
for i=1:a1
    if m(i)==rl
        m(i)=k;
    end
end
for i=1:a2
    if mm(i)==k
        mm(i)=rl;
    end
end
B=AA1(:,m);
%建立初始单纯形表；
cb=c1(m(1));
xb=(inv(B))*b;
for i=2:a1
    cb=[cb c1(m(i))];%取基向量cb;
end
z0=cb*xb;          %计算目标函数,即基本可行解；
rr=zeros(1,a1+a2);
%计算检验数；
for j=1:a2                          %计算检验数；
    cy=0;
    for i=1:a1
        cy=cy+c1(m(i))*AA1(i,mm(j));
    end
    z(mm(j))=cy;
    rr(mm(j))=c1(mm(j))-z(mm(j));
end
AA2=[rr -z0];
disp('各步单纯性表')
AA=[AA1;AA2]    %生成初始单纯形表；
if (min(rr))>=0    %判断是否是最优解，如果是则显示变量取值x以及最优值z；
    z=z0;
    x(m)=xb;
    disp('变量最优解为：x*=')
    disp(x);
    disp('最优解值为：z*=')
    disp(z);
    flag=0;
    break;
end
%决定进基矢量ak;
k=min(find(rr==min(rr(find(rr<0)))));     
ak=AA1(:,k);
%决定离基矢量ar和主元素yrk;
if max(AA1(:,k))<=0        %如果yik中有小于等于0的值则目标函数是无界的，无最优解；
    disp('无有限最优解')
    flag=0;
    break;
end  
for i=1:a1
    NN(i)=(AA1(i,a1+a2+1))/(AA1(i,k));
end
r=find(NN==min(NN(find(NN>0))));
    yrk=AA1(r,k);            %主元素为yrk； 
    rl=m(r);                %离基矢量下标为rl；
    ar=AA1(:,rl);         %离基矢量为ar；
%根据列主元Gauss消元法解线性方程组
for i=1:a1
    if i~=r
        AA1(i,:)=AA1(i,:)-AA1(r,:)*AA1(i,k)/AA1(r,k);
    else
        AA1(r,:)=AA1(r,:)/yrk;
    end
end
AA1;
end



   
    
    