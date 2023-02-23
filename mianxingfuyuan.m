clear all;
m=128;
n=160;

%***构造一组参考镜与参考镜干涉的四步移相干涉图***%
ff=n;
t=(0:ff-1)/ff;
A=1;
f=4;
coh1=2*A^2+2*A^2*cos(2*pi*f*t+0);
coh2=2*A^2+2*A^2*cos(2*pi*f*t+1.57);
coh3=2*A^2+2*A^2*cos(2*pi*f*t+3.14);
coh4=2*A^2+2*A^2*cos(2*pi*f*t+4.71);
coh5=2*A^2+2*A^2*cos(2*pi*f*t+6);
%注意，ccd5不可以取6.28，否则和ccd1差距太小，导致计算结果相对来说出现较大的误差
% 5次移相不要求全部递增，但是不要只局限在0-pi之间，至少要有一个位于pi-2pi之间，且差别越大，最终结果越精确。否则的话，误差可能会比较大。

%***设定初始值***%
delta=[0,1,2.5,4,5.5];%假设初始的相位移动量,单位弧度
%初始值要相对来说接近真值，否则最终结果也可能会出现较大偏差

%***构造CCD拍摄图***%
ccd1=zeros(m,n);
ccd2=zeros(m,n);
ccd3=zeros(m,n);
ccd4=zeros(m,n);
ccd5=zeros(m,n);
for mm=1:m
    ccd1(mm,:)=coh1;%ccd拍摄图
    ccd2(mm,:)=coh2;
    ccd3(mm,:)=coh3;
    ccd4(mm,:)=coh4;
    ccd5(mm,:)=coh5;
end;


%***初始化***%
J=size(delta,2);%步数
deltak=zeros(1,J);%用来存储上一个迭代算得的相位移动量,单位弧度
deltaD=abs(delta-deltak);%迭代之后的相位变化量,单位弧度
I=zeros(m,n,J);%光强图初始化
psi=I;%psi矩阵初始化
% mesh(ccd1);%输出三维图

%***构造光强矩阵***%
I(:,:,1)=ccd1;
I(:,:,2)=ccd2;
I(:,:,3)=ccd3;
I(:,:,4)=ccd4;
I(:,:,5)=ccd5;

%***迭代相关数值初始化***%
S=zeros(m,n);
C=zeros(m,n);
sindelta=zeros(1,J);
cosdelta=zeros(1,J);
tandelta=zeros(1,J);
lamda=570;%单位nm
phi=zeros(m,n);%面形所对应的相位
% fs=phi*lamda/(2*pi);%即face shape面形
kappa=1*10^(-3);
num=0;%迭代次数

%***构造psi矩阵***%
% psi(:,:,1)=I(:,:,1);
for j=1:J
    psi(:,:,j)=I(:,:,j)-I(:,:,1);%psi矩阵赋值
end;

%***迭代循环***%
while deltaD(1)>=kappa||deltaD(2)>=kappa||deltaD(3)>=kappa||deltaD(4)>=kappa||deltaD(5)>=kappa
    %***初始化***%
    a=0;b=0;c=0;%a,b,c赋初值
    f=0;g=0;h=0;%f,g,h赋初值
    d=zeros(m,n);
    e=zeros(m,n);
    s=zeros(1,J);
    t=zeros(1,J);
    
    %***连续最小二乘拟合计算S和C***%
    for j=1:J
        a=a+(cos(delta(j))-1)^2;
        b=b+sin(delta(j))*(cos(delta(j))-1);
        c=c+(sin(delta(j)))^2;
        for p=1:m
            for q=1:n
                d(p,q)=d(p,q)+psi(p,q,j)*(cos(delta(j))-1);
                e(p,q)=e(p,q)+psi(p,q,j)*sin(delta(j));
            end;
        end;
    end;
    for p=1:m
        for q=1:n
            S(p,q)=(a*e(p,q)-b*d(p,q))/(a*c-b^2);
            C(p,q)=(c*d(p,q)-b*e(p,q))/(a*c-b^2);
            phi(p,q)=atan(S(p,q)/C(p,q));
        end;
    end;
    
    %***空间最小二乘拟合计算delta***%
    for p=1:m
        for q=1:n
            f=f+(C(p,q))^2;
            g=g+C(p,q)*S(p,q);
            h=h+(S(p,q))^2;
            for j=1:J
                s(j)=s(j)+psi(p,q,j)*C(p,q)+(C(p,q))^2;
                t(j)=t(j)+psi(p,q,j)*S(p,q)+S(p,q)*C(p,q);
            end;
        end;
    end;
    deltak=delta;
    for j=1:J
        sindelta(j)=(f*t(j)-g*s(j))/(f*h-g^2);
        cosdelta(j)=(h*s(j)-g*t(j))/(f*h-g^2);
        if (sindelta(j)>=0&&cosdelta(j)>0)
            delta(j)=atan(sindelta(j)/cosdelta(j));
        elseif (sindelta(j)>=0&&cosdelta(j)<0)
            delta(j)=atan(sindelta(j)/cosdelta(j))+pi;
        elseif (sindelta(j)<=0&&cosdelta(j)<0)
            delta(j)=atan(sindelta(j)/cosdelta(j))+pi;
        elseif (sindelta(j)<0&&cosdelta(j)>0)
            delta(j)=atan(sindelta(j)/cosdelta(j))+2*pi;
        elseif (sindelta(j)>0&&cosdelta(j)==0)
            delta(j)=pi/2;
        elseif (sindelta(j)<0&&cosdelta(j)==0)
            delta(j)=3*pi/2;
        end;
    end;
    deltaD=abs(delta-deltak);
    num=num+1;
end;

%***计算最终的phi***%
for j=1:J
    a=a+(cos(delta(j))-1)^2;
    b=b+sin(delta(j))*(cos(delta(j))-1);
    c=c+(sin(delta(j)))^2;
    for p=1:m
        for q=1:n
            d(p,q)=d(p,q)+psi(p,q,j)*(cos(delta(j))-1);
            e(p,q)=e(p,q)+psi(p,q,j)*(sin(delta(j)));
        end;
    end;
end;
for p=1:m
    for q=1:n
        S(p,q)=(a*e(p,q)-b*d(p,q))/(a*c-b^2);
        C(p,q)=(c*d(p,q)-b*e(p,q))/(a*c-b^2);
        if (S(p,q)>=0&&C(p,q)>0)
            phi(p,q)=atan(S(p,q)/C(p,q));
        elseif (S(p,q)>=0&&C(p,q)<0)
            phi(p,q)=atan(S(p,q)/C(p,q))+pi;
        elseif (S(p,q)<=0&&C(p,q)<0)
            phi(p,q)=atan(S(p,q)/C(p,q))+pi;
        elseif (S(p,q)<0&&C(p,q)>0)
            phi(p,q)=atan(S(p,q)/C(p,q))+2*pi;
        elseif (S(p,q)>0&&C(p,q)==0)
            phi(p,q)=pi/2;
        elseif (S(p,q)<0&&C(p,q)==0)
            phi(p,q)=3*pi/2;
        end;
    end;
end;

fs=phi*lamda/(2*pi);
mesh(phi);%输出三维图