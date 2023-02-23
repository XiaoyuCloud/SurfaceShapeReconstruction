%较好数据：19、31、32、33、41、42、43
clc;
%***读取CCD拍摄图***%
J=9;jn=1;
xun=[1,2,3,4,5,6,7,8,9];
m1=494;
m2=648;
n1=468;
n2=612;
m=m2-m1+1;
n=n2-n1+1;
I=zeros(m,n,J);%光强图初始化
kappa=1*10^(-3);%即kappa
quyu=31;
for quyu=quyu:quyu
    for j=1:9
        pathname=strcat('.\11.18\',num2str(quyu),'\',num2str(j),'.jpg');
        ccdgray=rgb2gray(imread(pathname));
        I(:,:,jn)=ccdgray(m1:m2,n1:n2);
        jn=jn+1;
    end;
end;

%***设定初始值***%
delta=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
%假设初始的相位移动量,单位弧度
%初始值要相对来说接近真值，否则最终结果也可能会出现较大偏差
% delta=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5];
deltali=delta;

%***初始化***%
deltak=zeros(1,J);%用来存储上一个迭代算得的相位移动量,单位弧度
deltaD=abs(delta-deltak);%迭代之后的相位变化量,单位弧度
% I=zeros(m,n,J);%光强图初始化
psi=zeros(m,n,J);%psi矩阵初始化

%***迭代相关数值初始化***%
S=zeros(m,n);
C=zeros(m,n);
sindelta=zeros(1,J);
cosdelta=zeros(1,J);
tandelta=zeros(1,J);
lamda=570;%单位nm
phi=zeros(m,n);%面形所对应的相位
phichu=zeros(m,n);%atan直接得到的phi
num=0;%迭代次数

%***四步移相相关数据初始化***%
phifour=zeros(m,n);
sinfour=zeros(m,n);
cosfour=zeros(m,n);
deltafour=zeros(1,J);

%***构造psi矩阵***%
for j=1:J
    psi(:,:,j)=I(:,:,j)-I(:,:,1);%psi矩阵赋值
end;

%||deltaD(6)>=kappa||deltaD(7)>=kappa||deltaD(8)>=kappa||deltaD(9)>=kappa
%***迭代循环***%
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

% %***将phi从0变到2pi时的跳变去除掉***%
% for p=1:m
%     Q=1;
%     while Q<=(n-2)&&(phi(p,Q)<5||phi(p,Q+1)<5||phi(p,Q+2)<5)
%         Q=Q+1;
%     end;
%     for q=Q:-1:2
%         if((phi(p,q)-phi(p,q-1))>=5&&abs(I(p,q,1)-I(p,q-1,1))<=200)
%             phi(p,q-1)=phi(p,q-1)+2*pi;
%         end;
%     end;
%     for q=Q:n-1
%         if((phi(p,q)-phi(p,q+1))>=5&&abs(I(p,q,1)-I(p,q+1,1))<=200)
%             phi(p,q+1)=phi(p,q+1)+2*pi;
%         end;
%     end;
% end;
% for q=1:n
%     P=1;
%     while P<=(m-2)&&(phi(P,q)<6||phi(P+1,q)<6||phi(P+2,q)<6)
%         P=P+1;
%     end;
%     for p=P:-1:2
%         if((phi(p,q)-phi(p-1,q))>=5&&abs(I(p,q,1)-I(p-1,q,1))<=200)
%             phi(p-1,q)=phi(p-1,q)+2*pi;
%         end;
%     end;
%     for p=P:m-1
%         if((phi(p,q)-phi(p+1,q))>=5&&abs(I(p,q,1)-I(p+1,q,1))<=200)
%             phi(p+1,q)=phi(p+1,q)+2*pi;
%         end;
%     end;
% end;

fs=0.5*phi*lamda/(2*pi);
% axis([1 m 1 n 400 500]);
deltabi=[deltali;delta;deltali-delta];

% figure(1);
% imshow(ccd5part);
% x=1:9;
% dian=[131,87];
% white=[ccd1part(dian(1,1),dian(1,2)),ccd2part(dian(1,1),dian(1,2)),ccd3part(dian(1,1),dian(1,2)),ccd4part(dian(1,1),dian(1,2)),ccd5part(dian(1,1),dian(1,2)),ccd6part(dian(1,1),dian(1,2)),ccd7part(dian(1,1),dian(1,2)),ccd8part(dian(1,1),dian(1,2)),ccd9part(dian(1,1),dian(1,2))];
% figure(2);
% plot(x,white);

%***四步移相***%
deltaza=[0,0.7,1.4,2.1,2.8,3.5,4.2,4.9,5.6];
deltafour=deltali;%结论，根初始的deltali或者最小二乘拟合得到的delta，用四步移向法计算和用最小二乘法计算的面形之间,相关系数在0.999以上才比较好。如果数据相差较大，相关系数并不会下降很多，所以从相关系数并不能看出来什么
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

% %***将phi从0变到2pi时的跳变去除掉***%
% for p=1:m
%     Q=1;
%     while Q<=(n-2)&&(phifour(p,Q)<5||phifour(p,Q+1)<5||phifour(p,Q+2)<5)
%         Q=Q+1;
%     end;
%     for q=Q:-1:2
%         if((phifour(p,q)-phifour(p,q-1))>=5&&abs(I(p,q,1)-I(p,q-1,1))<=200)
%             phifour(p,q-1)=phifour(p,q-1)+2*pi;
%         end;
%     end;
%     for q=Q:n-1
%         if((phifour(p,q)-phifour(p,q+1))>=5&&abs(I(p,q,1)-I(p,q+1,1))<=200)
%             phifour(p,q+1)=phifour(p,q+1)+2*pi;
%         end;
%     end;
% end;
% for q=1:n
%     P=1;
%     while P<=(m-2)&&(phifour(P,q)<6||phifour(P+1,q)<6||phifour(P+2,q)<6)
%         P=P+1;
%     end;
%     for p=P:-1:2
%         if((phifour(p,q)-phifour(p-1,q))>=5&&abs(I(p,q,1)-I(p-1,q,1))<=200)
%             phifour(p-1,q)=phifour(p-1,q)+2*pi;
%         end;
%     end;
%     for p=P:m-1
%         if((phifour(p,q)-phifour(p+1,q))>=5&&abs(I(p,q,1)-I(p+1,q,1))<=200)
%             phifour(p+1,q)=phifour(p+1,q)+2*pi;
%         end;
%     end;
% end;

figure(1)
subplot(1,2,1);
surf(phi);%输出三维图
axis([1 m 1 n -pi 2*pi]);
shading interp;
subplot(1,2,2);
surf(phifour);%输出三维图
shading interp;
axis([1 m 1 n -pi 2*pi]);

Rdelta=corrcoef(delta,deltafour);
Rphi=corrcoef(phi,phifour);