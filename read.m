clc;

%***��ȡCCD����ͼ***%
ccd1=imread('.\11.18\12\1.jpg');
ccd2=imread('.\11.18\12\2.jpg');
ccd3=imread('.\11.18\12\3.jpg');
ccd4=imread('.\11.18\12\4.jpg');
ccd5=imread('.\11.18\12\5.jpg');

%***��CCD����ͼתΪ�Ҷ�ͼ***%
ccd1gray=rgb2gray(ccd1);
ccd2gray=rgb2gray(ccd2);
ccd3gray=rgb2gray(ccd3);
ccd4gray=rgb2gray(ccd4);
ccd5gray=rgb2gray(ccd5);
% imshow(ccd5gray);

%***������ͼ֮������򶼱�Ϊ0***%
% ccd1gray(1:470,n)=0;ccd2gray(1:470,n)=0;ccd3gray(1:470,n)=0;ccd4gray(1:470,n)=0;ccd5gray(1:470,n)=0;
% ccd1gray(780:m,n)=0;ccd2gray(780:m,n)=0;ccd3gray(780:m,n)=0;ccd4gray(780:m,n)=0;ccd5gray(780:m,n)=0;
% ccd1gray(470:780,1:490)=0;ccd2gray(470:780,1:490)=0;ccd3gray(470:780,1:490)=0;ccd4gray(470:780,1:490)=0;ccd5gray(470:780,1:490)=0;
% ccd1gray(470:780,850:n)=0;ccd2gray(470:780,850:n)=0;ccd3gray(470:780,850:n)=0;ccd4gray(470:780,850:n)=0;ccd5gray(470:780,850:n)=0;
% imshow(ccd1gray);
% ccd1part=ccd1gray(355:645,505:795);
% ccd2part=ccd2gray(355:645,505:795);
% ccd3part=ccd3gray(355:645,505:795);
% ccd4part=ccd4gray(355:645,505:795);
% ccd5part=ccd5gray(355:645,505:795);

ccd1part=ccd1gray(400:600,550:750);
ccd2part=ccd2gray(400:600,550:750);
ccd3part=ccd3gray(400:600,550:750);
ccd4part=ccd4gray(400:600,550:750);
ccd5part=ccd5gray(400:600,550:750);
% imshow(ccd1part);
m=size(ccd1part,1);
n=size(ccd2part,2);

%***�趨��ʼֵ***%
delta=[0,2*pi/5,4*pi/5,6*pi/5,8*pi/5];%�����ʼ����λ�ƶ���,��λ����
%��ʼֵҪ�����˵�ӽ���ֵ���������ս��Ҳ���ܻ���ֽϴ�ƫ��

%***��ʼ��***%
J=size(delta,2);%����
deltak=zeros(1,J);%�����洢��һ��������õ���λ�ƶ���,��λ����
deltaD=abs(delta-deltak);%����֮�����λ�仯��,��λ����
I=zeros(m,n,J);%��ǿͼ��ʼ��
psi=I;%psi�����ʼ��
% mesh(ccd1);%�����άͼ

%***��ǿ����ֵ***%
I(:,:,1)=ccd1part;
I(:,:,2)=ccd2part;
I(:,:,3)=ccd3part;
I(:,:,4)=ccd4part;
I(:,:,5)=ccd5part;

%***���������ֵ��ʼ��***%
S=zeros(m,n);
C=zeros(m,n);
sindelta=zeros(1,J);
cosdelta=zeros(1,J);
tandelta=zeros(1,J);
lamda=570;%��λnm
phi=zeros(m,n);%��������Ӧ����λ
% fs=phi*lamda/(2*pi);%��face shape����
kappa=1*10^(-2);
num=0;%��������

%***����psi����***%
% psi(:,:,1)=I(:,:,1);
for j=1:J
    psi(:,:,j)=I(:,:,j)-I(:,:,1);%psi����ֵ
end;

%***����ѭ��***%
while deltaD(1)>=kappa||deltaD(2)>=kappa||deltaD(3)>=kappa||deltaD(4)>=kappa||deltaD(5)>=kappa
    %***��ʼ��***%
    a=0;b=0;c=0;%a,b,c����ֵ
    f=0;g=0;h=0;%f,g,h����ֵ
    d=zeros(m,n);
    e=zeros(m,n);
    s=zeros(1,J);
    t=zeros(1,J);
    
    %***������С������ϼ���S��C***%
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
    
    %***�ռ���С������ϼ���delta***%
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

%***�������յ�phi***%
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

%***��phi��0�䵽2piʱ������ȥ����***%
for p=1:m
    Q=1;
    while Q<=(n-2)&&(phi(p,Q)<5||phi(p,Q+1)<5||phi(p,Q+2)<5)
        Q=Q+1;
    end;
    for q=Q:-1:2
        if((phi(p,q)-phi(p,q-1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q-1))<=300)
            phi(p,q-1)=phi(p,q-1)+2*pi;
        end;
    end;
    for q=Q:n-1
        if((phi(p,q)-phi(p,q+1))>=5&&abs(ccd1part(p,q)-ccd1part(p,q+1))<=300)
            phi(p,q+1)=phi(p,q+1)+2*pi;
        end;
    end;
end;
for q=1:n
    P=1;
    while P<=(m-2)&&(phi(P,q)<6||phi(P+1,q)<6||phi(P+2,q)<6)
        P=P+1;
    end;
    for p=P:-1:2
        if((phi(p,q)-phi(p-1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p-1,q))<=300)
            phi(p-1,q)=phi(p-1,q)+2*pi;
        end;
    end;
    for p=P:m-1
        if((phi(p,q)-phi(p+1,q))>=5&&abs(ccd1part(p,q)-ccd1part(p+1,q))<=300)
            phi(p+1,q)=phi(p+1,q)+2*pi;
        end;
    end;
end;

fs=phi*lamda/(2*pi);
surf(fs);%�����άͼ
shading interp;
% axis([1 m 1 n 400 500]);