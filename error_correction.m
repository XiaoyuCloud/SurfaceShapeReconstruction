% 较好数据：19、41、42
clear all;
clc;

%***设定一些初始值***%
J=9;
xun=[1,2,3,4,5,6,7,8,9];%选择哪几幅干涉图进行计算
m1=1;
m2=200;
n1=1;
n2=200;
m=m2-m1+1;
n=n2-n1+1;
mn=110;
mcentre=round(m/2);
ncentre=round(n/2);
I=zeros(m,n,J);%I为干涉图所对应的灰度图，这里将其初始化
kappa=1*10^(-4);%即kappa，收敛判据
quyu1=81;%选择所测量的面形区域。XY意为第X块区域第Y次测量所得到的数据
quyu2=81;
Sumnum=quyu2-quyu1+1;
CharacterSum=zeros(Sumnum*5,9);%5*3
number=1;
lamda=570;%单位nm
phiquyu=zeros(m,n);%存储同一区域不同组数据减去零点并平均后的结果;




refer=1;
meaed=1;
psimax=0.5*pi;
hphimax=20*pi;
hdelta=hphimax/n;
R=0.5*psimax+mcentre^2/(2*psimax);
%***构造一个球面***%
Hphi=zeros(m,n);
phili=zeros(m,n);
for q=1:m
%     for p=1:n
%         Hphi(p,q)=R-(R^2-((p-mcentre)^2+(q-ncentre)^2))^(0.5);
%     end;
    Hphi(:,q)=Hphi(:,q)+(q-1)*hdelta;
end;
H=0.5*Hphi/(2*pi)*lamda;
deltaHphi=pi/5;
deltaH=0.5*deltaHphi/(2*pi)*lamda;%每次参考相位移动的deltah；

%***读取拍摄的干涉图***%
for quyu=quyu1:quyu2
    jn=1;%J为移相的次数。
    ccdnum=zeros(m,n,J);
    cohnum=zeros(J,n);
    for j=1:9
        for p=1:m
            for q=1:n
%                 ccdnum(p,q,j)=abs(meaed+refer*cos(2*(H(p,q)-(j-1)*deltaH)*2*pi/lamda));%9次移相  
%                 phili(p,q)=2*H(p,q)*2*pi/lamda;
                ccdnum(p,q,j)=abs(meaed+refer*cos(Hphi(p,q)-(j-1)*deltaHphi));%9次移相  
                phili(p,q)=Hphi(p,q);
            end;
        end;
    end;
    for j=1:9
        I(:,:,jn)=ccdnum(:,:,jn);
        jn=jn+1;
    end;
    ccd1part=I(:,:,1);
    figure(15)
    for j=1:9
        imshow(I(:,:,j))
    end;

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
    phi=zeros(m,n);%被测面形所对应的相位
    phichu=zeros(m,n);%atan直接得到的phi
    phih=zeros(m,n);
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
                phi(p,q)=atan(S(p,q)/C(p,q));
            elseif (S(p,q)>0&&C(p,q)==0)
                phi(p,q)=pi/2;
            elseif (S(p,q)<0&&C(p,q)==0)
                phi(p,q)=3*pi/2;
            end;
        end;
    end;
    %***相位解包***%
    for numsin=1:8 %与最大相位差有关
        %mcentre那一行横向
        p=mcentre;
        K=0;
        qsingul=zeros(1,10);
        cb=0.5;
        for q=ncentre:-1:2
            if(abs(phi(p,q)-phi(p,q-1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q-1))<=cb)
                K=K+1;
                qsingul(K)=q-1;
                sym=(phi(p,q)-phi(p,q-1))/abs(phi(p,q)-phi(p,q-1));
            end;
        end;
        if(K~=0)
            for k=1:K
                if(k~=K)
                    phi(p,(qsingul(k+1)+1):qsingul(k))=phi(p,(qsingul(k+1)+1):qsingul(k))+k*2*pi*sym;
                else
                    phi(p,1:qsingul(k))=phi(p,1:qsingul(k))+k*2*pi*sym;
                end;
            end;
        end;

        K=0;
        qsingul=zeros(1,10);
        for q=ncentre:1:n-1
            if(abs(phi(p,q)-phi(p,q+1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q+1))<=cb)
                K=K+1;
                qsingul(K)=q+1;
                sym=(phi(p,q)-phi(p,q+1))/abs(phi(p,q)-phi(p,q+1));
            end;
        end;
        if(K~=0)
            for k=1:K
                if(k~=K)
                    phi(p,qsingul(k):(qsingul(k+1)-1))=phi(p,qsingul(k):(qsingul(k+1)-1))+k*2*pi*sym;
                else
                    phi(p,qsingul(k):n)=phi(p,qsingul(k):n)+k*2*pi*sym;
                end;
            end;
        end;

        %纵向
        for q=1:n
            K=0;
            psingul=zeros(1,10);
            sym=0;
            for p=mcentre:-1:2
                if(abs(phi(p,q)-phi(p-1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p-1,q))<=cb)
                    K=K+1;
                    psingul(K)=p-1;
                    sym=(phi(p,q)-phi(p-1,q))/abs(phi(p,q)-phi(p-1,q));
                end;
            end;
            if(K~=0)
                for k=1:K
                    if(k~=K)
                        phi((psingul(k+1)+1):psingul(k),q)=phi((psingul(k+1)+1):psingul(k),q)+k*2*pi*sym;
                    else
                        phi(1:psingul(k),q)=phi(1:psingul(k),q)+k*2*pi*sym;
                    end;
                end;
            end;

            K=0;
            psingul=zeros(1,10);
            sym=0;
            for p=mcentre:1:m-1
                if(abs(phi(p,q)-phi(p+1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p+1,q))<=cb)
                    K=K+1;
                    psingul(K)=p+1;
                    sym=(phi(p,q)-phi(p+1,q))/abs(phi(p,q)-phi(p+1,q));
                end;
            end;
            if(K~=0)
                for k=1:K
                    if(k~=K)
                        phi(psingul(k):(psingul(k+1)-1),q)=phi(psingul(k):(psingul(k+1)-1),q)+k*2*pi*sym;
                    else
                        phi(psingul(k):m,q)=phi(psingul(k):m,q)+k*2*pi*sym;
                    end;
                end;
            end;
        end;
    end;
    
    %***减去偏置***%
    fs=0.5*phi*lamda/(2*pi);%即face shape
    deltabi=[deltali;deltachu;delta];
    
    X=[];
    for i=1:n
         iI=i*ones(1,n);
         X=[X;iI]; 
    end
    X=X';
    X=X(:);
    Y=[];
    for i=1:m
        iI=i*ones(1,m);
        Y=[Y;iI];
    end
    Y=Y(:);
    Z=phi(:);
    phipoly=fit([X,Y],Z,'poly11');
    phipolyZ=phipoly(X,Y)';
    phiplane=zeros(m,n);
    phiplanexie=zeros(m,n);
    tidu=(max(phi(:))-min(phi(:)))/(n-1);
    for p=1:m
            phiplane(:,p)=(phi(:,p)-(phipolyZ(1,((p-1)*m+1):p*m))');
            phiplanexie(:,p)=(phi(:,p)-(p-1)*tidu);
    end;
    
    phiquyumin=min(phi(:));
    phiquyu=phiquyu+phi-phiquyumin;
    
    %***四步移相算法***%
    deltafour=deltali;%结论，根初始的deltali或者最小二乘拟合得到的delta，用四步移向法计算和用最小二乘法计算的面形之间,相关系数在0.999以上才比较好。如果数据相差较大，相关系数并不会下降很多，所以从相关系数并不能看出来什么
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
                    phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q));
                elseif (sinfour(p,q)>0&&cosfour(p,q)==0)
                    phifour(p,q)=pi/2;
                elseif (sinfour(p,q)<0&&cosfour(p,q)==0)
                    phifour(p,q)=3*pi/2;
                end;
            end;
        end;
    end;

    phizanshi=phi;%暂时存储一下phi的数据
    phi=phifour;%下面代码就不用修改了
    %***九步移相相位解包***%
    %mcentre那一行横向
    for ibi=1:8
        p=mcentre;
        K=0;
        qsingul=zeros(1,10);
        cb=100;
        for q=ncentre:-1:2
            if(abs(phi(p,q)-phi(p,q-1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q-1))<=cb)
                K=K+1;
                qsingul(K)=q-1;
                sym=(phi(p,q)-phi(p,q-1))/abs(phi(p,q)-phi(p,q-1));
            end;
        end;
        if(K~=0)
            for k=1:K
                if(k~=K)
                    phi(p,(qsingul(k+1)+1):qsingul(k))=phi(p,(qsingul(k+1)+1):qsingul(k))+k*2*pi*sym;
                else
                    phi(p,1:qsingul(k))=phi(p,1:qsingul(k))+k*2*pi*sym;
                end;
            end;
        end;

        K=0;
        qsingul=zeros(1,10);
        for q=ncentre:1:n-1
            if(abs(phi(p,q)-phi(p,q+1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q+1))<=cb)
                K=K+1;
                qsingul(K)=q+1;
                sym=(phi(p,q)-phi(p,q+1))/abs(phi(p,q)-phi(p,q+1));
            end;
        end;
        if(K~=0)
            for k=1:K
                if(k~=K)
                    phi(p,qsingul(k):(qsingul(k+1)-1))=phi(p,qsingul(k):(qsingul(k+1)-1))+k*2*pi*sym;
                else
                    phi(p,qsingul(k):n)=phi(p,qsingul(k):n)+k*2*pi*sym;
                end;
            end;
        end;

        %纵向
        for q=1:n
            K=0;
            psingul=zeros(1,10);
            sym=0;
            for p=mcentre:-1:2
                if(abs(phi(p,q)-phi(p-1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p-1,q))<=cb)
                    K=K+1;
                    psingul(K)=p-1;
                    sym=(phi(p,q)-phi(p-1,q))/abs(phi(p,q)-phi(p-1,q));
                end;
            end;
            if(K~=0)
                for k=1:K
                    if(k~=K)
                        phi((psingul(k+1)+1):psingul(k),q)=phi((psingul(k+1)+1):psingul(k),q)+k*2*pi*sym;
                    else
                        phi(1:psingul(k),q)=phi(1:psingul(k),q)+k*2*pi*sym;
                    end;
                end;
            end;

            K=0;
            psingul=zeros(1,10);
            sym=0;
            for p=mcentre:1:m-1
                if(abs(phi(p,q)-phi(p+1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p+1,q))<=cb)
                    K=K+1;
                    psingul(K)=p+1;
                    sym=(phi(p,q)-phi(p+1,q))/abs(phi(p,q)-phi(p+1,q));
                end;
            end;
            if(K~=0)
                for k=1:K
                    if(k~=K)
                        phi(psingul(k):(psingul(k+1)-1),q)=phi(psingul(k):(psingul(k+1)-1),q)+k*2*pi*sym;
                    else
                        phi(psingul(k):m,q)=phi(psingul(k):m,q)+k*2*pi*sym;
                    end;
                end;
            end;
        end;
%         if ibi==1
%             phiweiqu1=phi;
%         end;
    end;
    phifour=phi;
    phi=phizanshi;

    picphi=figure(1);
    contourf(phi);%输出三维图
    title('最小二乘拟合迭代算法');
%     axis([1 m 1 n 0 6*pi]);

    picphiplane=figure(2);
    contourf(phiplane);
%     axis([1 m 1 n -1 2]); 
    title('最小二乘拟合-减去等倾干涉偏置');
    
    picphiplanexie=figure(3);
    contourf(phiplanexie);
%   axis([1 m 1 n -1 2]); 
    title('最小二乘拟合-减去斜面等倾干涉偏置');
    
    picphifour=figure(4);
    contourf(phifour);%输出三维图
%     axis([1 m 1 n 0 6*pi]);
    title('九步移相算法');
    
    figure(5);
    mesh(phi);%输出三维图
    title('最小二乘拟合迭代算法');
%     axis([1 m 1 n 0 6*pi]);

    figure(6);
    mesh(phiplane);
%     axis([1 m 1 n -1 2]); 
    title('最小二乘拟合-减去等倾干涉偏置');
        
    figure(7);
    mesh(phiplanexie);
%   axis([1 m 1 n -1 2]); 
    title('最小二乘拟合-减去斜面等倾干涉偏置');
    
    figure(8);
    mesh(phifour);%输出三维图
%     axis([1 m 1 n 0 6*pi]);
    title('九步移相算法');
    
%     
%     
%     figure(6);
%     mesh(phiweiqu1);
    
    namephi=strcat('phi',num2str(quyu),'.jpg');
    namephiplane=strcat('phiplane',num2str(quyu),'.jpg');
    namematphiplane=strcat(num2str(quyu),'.mat');
    namematphiplanexie=strcat('xie',num2str(quyu),'.mat');
    namephifour=strcat('phifour',num2str(quyu),'.jpg');
    saveas(picphi,namephi);
    saveas(picphiplane,namephiplane);
    save(namematphiplane,'phiplane');
    save(namematphiplanexie,'phiplanexie');
    saveas(picphifour,namephifour);
    
    phi=phiplane;
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
    CharacterSum(number:number+4,:)=[deltabi;character];
    number=number+5;
    
    figure(14)
    for j=1:9
        subplot(2,5,j);
        ccd1part=I(:,:,j);%方便下面代码不用修改了
        imshow(ccd1part);
    end;
end;
nameCS=strcat(num2str(quyu1),'_',num2str(quyu2),'.xlsx');
xlswrite(nameCS,CharacterSum);

phiquyu=phiquyu/(quyu2-quyu1+1);
figure(9);
contourf(phiquyu);
title('最小二乘拟合平均并减去零点');
namematphiquyu=strcat('quyu',num2str(floor(quyu2/10)),'.mat');
save(namematphiquyu,'phiquyu');

figure(10);
mesh(phiquyu);
title('最小二乘拟合平均并减去零点');

X=[];
for i=1:n
     iI=i*ones(1,n);
     X=[X;iI]; 
end
X=X';
X=X(:);
Y=[];
for i=1:m
    iI=i*ones(1,m);
    Y=[Y;iI];
end
Y=Y(:);
Z=phiquyu(:);
phiquyupoly=fit([X,Y],Z,'poly11');
phiquyupolyZ=phiquyupoly(X,Y)';
phiquyuplane=zeros(m,n);
phiquyuplanexie=zeros(m,n);
tiduquyu=(max(phiquyu(:))-min(phiquyu(:)))/(n-1);
for p=1:m
        phiquyuplane(:,p)=(phiquyu(:,p)-(phiquyupolyZ(1,((p-1)*m+1):p*m))');
        phiquyuplanexie(:,p)=(phiquyu(:,p)-(p-1)*tiduquyu);
end;
% for q=1:n
%     phiquyuplane(:,q)=phiquyu(:,q)-(q-1)*hdelta;
%     phiquyuplanexie(:,p)=(phiquyu(:,p)-(m-p+1)*tiduquyu);
% end;
picphiquyuplane=figure(11);
mesh(phiquyuplane);
title('最小二乘拟合平均并减去零点减去拟合偏置');
quyu=floor(quyu2/10);
namephiquyuplane=strcat('phiquyuplane',num2str(quyu),'.jpg');
saveas(picphiquyuplane,namephiquyuplane);
figure(12);
mesh(phiquyuplanexie);
title('最小二乘拟合平均并减去零点减去直接偏置');
figure(13)
plot(phiquyupoly,[X,Y],Z);