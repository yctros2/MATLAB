function[A,N,n,a,b,m,flag,mA,nA]=simplex(A,m,b,c)%AΪ������NΪA�����һ�У�nΪ��ʼ����������aΪ������bΪ������mΪ����������flagΪcj��־��mA��nAΪA������������%�˴������С��ֵ��
% ����:
%   A - ��ʼϵ������
%   m - ������
%   b - Լ�������Ҳ�����
%   c - Ŀ�꺯��ϵ������
% ���:
%   A, N, n, a, b, m, flag, mA, nA - �����α���ز���
A=[3 -1 2;-2 4 0;-4 3 8];
b=[7;12;12];
c=[1 -3 -1];
A= sym(A);
b = sym(b);
c = sym(c);

 %���ɱ�׼�ͣ�
[a1 a2]=size(A);    %ȡA������a1,����a2��
B=eye(a1);        %ȡ��������ȵĵ�λ��
A1=[A B];         %��׼�������Ա���ϵ������
m=[(a2+1):(a1+a2)]; %�г����������±�ţ�
mm=[1:a2];         %�г��ǻ������±ꣻ
m1=zeros(1,((a1+a2)-size(c,2)));%��c����0������
c1=[c m1];                %����������c������
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
%������ʼ�����α�
cb=c1(m(1));
xb=(inv(B))*b;
for i=2:a1
    cb=[cb c1(m(i))];%ȡ������cb;
end
z0=cb*xb;          %����Ŀ�꺯��,���������н⣻
rr=zeros(1,a1+a2);
%�����������
for j=1:a2                          %�����������
    cy=0;
    for i=1:a1
        cy=cy+c1(m(i))*AA1(i,mm(j));
    end
    z(mm(j))=cy;
    rr(mm(j))=c1(mm(j))-z(mm(j));
end
AA2=[rr -z0];
disp('���������Ա�')
AA=[AA1;AA2]    %���ɳ�ʼ�����α�
if (min(rr))>=0    %�ж��Ƿ������Ž⣬���������ʾ����ȡֵx�Լ�����ֵz��
    z=z0;
    x(m)=xb;
    disp('�������Ž�Ϊ��x*=')
    disp(x);
    disp('���Ž�ֵΪ��z*=')
    disp(z);
    flag=0;
    break;
end
%��������ʸ��ak;
k=min(find(rr==min(rr(find(rr<0)))));     
ak=AA1(:,k);
%�������ʸ��ar����Ԫ��yrk;
if max(AA1(:,k))<=0        %���yik����С�ڵ���0��ֵ��Ŀ�꺯�����޽�ģ������Ž⣻
    disp('���������Ž�')
    flag=0;
    break;
end  
for i=1:a1
    NN(i)=(AA1(i,a1+a2+1))/(AA1(i,k));
end
r=find(NN==min(NN(find(NN>0))));
    yrk=AA1(r,k);            %��Ԫ��Ϊyrk�� 
    rl=m(r);                %���ʸ���±�Ϊrl��
    ar=AA1(:,rl);         %���ʸ��Ϊar��
%��������ԪGauss��Ԫ�������Է�����
for i=1:a1
    if i~=r
        AA1(i,:)=AA1(i,:)-AA1(r,:)*AA1(i,k)/AA1(r,k);
    else
        AA1(r,:)=AA1(r,:)/yrk;
    end
end
AA1;
end



   
    
    