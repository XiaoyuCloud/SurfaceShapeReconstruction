%较好数据：19、41、42
clc;

%***设定一些初始值***%
J=9;
xun=[1,2,3,4,5,6,7,8,9];%选择哪几幅干涉图进行计算
m1=494;%m1、m2、m3、m4用来选取相机所拍摄的图像中干涉图所在的区域。
m2=648;
n1=468;
n2=612;
m=m2-m1+1;
n=n2-n1+1;
I=zeros(m,n,J);%I为干涉图所对应的灰度图，这里将其初始化
kappa=1*10^(-3);%即kappa，收敛判据
quyu1=81;%选择所测量的面形区域。XY意为第X块区域第Y次测量所得到的数据
quyu2=81;
Sumnum=quyu2-quyu1+1;
CharacterSum=zeros(Sumnum*6,9);%5*3
number=1;

ff=n;
r=(0:ff-1)/ff;
A=1;
f=1;

%***读取拍摄的干涉图***%
for quyu=quyu1:quyu2
    jn=1;%J为移相的次数。
    ccdnum=zeros(m,n,J);
    cohnum=zeros(J,n);
    for j=1:9
        cohnum(j,:)=2*A^2+2*A^2*cos(2*pi*f*r+(j-1)*pi/5);
    end;
%     for theta=(0:2*pi)/ff
%         ccdnum(round())
%     end;
    for j=1:9
        for mm=1:m
            ccdnum(mm,:,j)=cohnum(j,:);%ccd拍摄图
        end;
    end;
    for j=1:9
        I(:,:,jn)=ccdnum(:,:,jn);
        jn=jn+1;
    end;
%     imshow(I(:,:,3));

    %***设定参考相位移动的初始值***%
%     delta=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
    delta=[0,0.1,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5];
      %假设初始的相位移动量,单位弧度
      %初始值要相对来说接近真值，否则最终结果也可能会出现较大偏差
    deltachu=delta;
    deltali=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
    
    %***初始化***%
    deltak=zeros(1,J);%用来存储上一个迭代算得的相位移动量,单位弧度
    deltaD=abs(delta-deltak);%迭代之后的相位变化量,单位弧度
    psi=zeros(m,n,J);%psi矩阵初始化

    %***迭代相关数值初始化***%
    S=zeros(m,n);
    C=zeros(m,n);
    sindelta=zeros(1,J);
    cosdelta=zeros(1,J);
    tandelta=zeros(1,J);
    lamda=570;%单位nm
    phi=zeros(m,n);%被测面形所对应的相位
    phichu=zeros(m,n);%atan直接得到的phi
    num=0;%迭代次数

    %***四步移相相关数据初始化***%
    phifour=zeros(m,n);%四步移相最后算得的面形所对应的相位
    sinfour=zeros(m,n);
    cosfour=zeros(m,n);

    %***构造psi矩阵***%
    for j=1:J
        psi(:,:,j)=I(:,:,j)-I(:,:,1);%psi矩阵赋值
    end;

    %***最小二乘拟合迭代算法***%
    while deltaD(1)>=kappa||deltaD(2)>=kappa||deltaD(3)>=kappa||deltaD(4)>=kappa||deltaD(5)>=kappa||deltaD(6)>=kappa||deltaD(7)>=kappa||deltaD(8)>=kappa||deltaD(9)>=kappa
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

    %***计算相应的phi***%
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
            phichu(p,q)=atan(S(p,q)/C(p,q));
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

    %***减去偏置***%
%     phimin=min(phi(:));
%     phi=phi-phimin;
    fs=0.5*phi*lamda/(2*pi);%即face shape
    deltabi=[deltali;deltachu;delta;deltali-delta];

    %***四步移相算法***%
    deltafour=deltachu;%结论，根初始的deltali或者最小二乘拟合得到的delta，用四步移向法计算和用最小二乘法计算的面形之间,相关系数在0.999以上才比较好。如果数据相差较大，相关系数并不会下降很多，所以从相关系数并不能看出来什么
                      %deltafour表示四步移相所用的参考相位
    for p=1:m
        for q=1:n
            for j=1:J
                sinfour(p,q)=sinfour(p,q)+I(p,q,j)*sin(deltafour(j));
                cosfour(p,q)=cosfour(p,q)+I(p,q,j)*cos(deltafour(j));
                if (sinfour(p,q)>=0&&cosfour(p,q)>0)
                   phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q));
                elseif (sinfour(p,q)>=0&&cosfour(p,q)<0)
                    phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q))+pi;
                elseif (sinfour(p,q)<=0&&cosfour(p,q)<0)
                    phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q))+pi;
                elseif (sinfour(p,q)<0&&cosfour(p,q)>0)
                    phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q))+2*pi;
                elseif (sinfour(p,q)>0&&cosfour(p,q)==0)
                    phifour(p,q)=pi/2;
                elseif (sinfour(p,q)<0&&cosfour(p,q)==0)
                    phifour(p,q)=3*pi/2;
                end;
            end;
        end;
    end;

    %***减去偏置***%
    phifourmin=min(phifour(:));
    phifour=phifour-phifourmin;

    picphi=figure(1);
    mesh(phi);%输出三维图
    title('最小二乘拟合迭代算法');
%     axis([1 m 1 n 0 2*pi]);
    picphifour=figure(2);
    mesh(phifour);%输出三维图
    title('九步移相算法');
    namephi=strcat('phi',num2str(quyu),'.jpg');
    namephifour=strcat('phifour',num2str(quyu),'.jpg');
    saveas(picphi,namephi);
    saveas(picphifour,namephifour);
    
    Rdelta=corrcoef(delta,deltafour);
    Rphi=corrcoef(phi,phifour);%迭代算法和移相算法最后所得的面形的相关系数
    MEANphi=mean(phi(:));
    MEANphifour=mean(phifour(:));
    RMSphi=rms(phi(:));
    RMSphifour=rms(phifour(:));
    MAXMINphi=max(phi(:))-min(phi(:));
    MAXMINphifour=max(phifour(:))-min(phifour(:));
    VARphi=var(phi(:));
    VARphifour=var(phifour(:));
    character=[MEANphi,MAXMINphi,RMSphi,VARphi,0,0,0,0,0;MEANphifour,MAXMINphifour,RMSphifour,VARphifour,0,0,0,0,0];
    CharacterSum(number:number+5,:)=[deltabi;character];
    number=number+6;
end;
nameCS=strcat(num2str(quyu1),'_',num2str(quyu2),'.xlsx');
xlswrite(nameCS,CharacterSum);
