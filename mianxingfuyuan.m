clear all;
m=128;
n=160;

%***����һ��ο�����ο���������Ĳ��������ͼ***%
ff=n;
t=(0:ff-1)/ff;
A=1;
f=4;
coh1=2*A^2+2*A^2*cos(2*pi*f*t+0);
coh2=2*A^2+2*A^2*cos(2*pi*f*t+1.57);
coh3=2*A^2+2*A^2*cos(2*pi*f*t+3.14);
coh4=2*A^2+2*A^2*cos(2*pi*f*t+4.71);
coh5=2*A^2+2*A^2*cos(2*pi*f*t+6);
%ע�⣬ccd5������ȡ6.28�������ccd1���̫С�����¼����������˵���ֽϴ�����
% 5�����಻Ҫ��ȫ�����������ǲ�Ҫֻ������0-pi֮�䣬����Ҫ��һ��λ��pi-2pi֮�䣬�Ҳ��Խ�����ս��Խ��ȷ������Ļ��������ܻ�Ƚϴ�

%***�趨��ʼֵ***%
delta=[0,1,2.5,4,5.5];%�����ʼ����λ�ƶ���,��λ����
%��ʼֵҪ�����˵�ӽ���ֵ���������ս��Ҳ���ܻ���ֽϴ�ƫ��

%***����CCD����ͼ***%
ccd1=zeros(m,n);
ccd2=zeros(m,n);
ccd3=zeros(m,n);
ccd4=zeros(m,n);
ccd5=zeros(m,n);
for mm=1:m
    ccd1(mm,:)=coh1;%ccd����ͼ
    ccd2(mm,:)=coh2;
    ccd3(mm,:)=coh3;
    ccd4(mm,:)=coh4;
    ccd5(mm,:)=coh5;
end;


%***��ʼ��***%
J=size(delta,2);%����
deltak=zeros(1,J);%�����洢��һ��������õ���λ�ƶ���,��λ����
deltaD=abs(delta-deltak);%����֮�����λ�仯��,��λ����
I=zeros(m,n,J);%��ǿͼ��ʼ��
psi=I;%psi�����ʼ��
% mesh(ccd1);%�����άͼ

%***�����ǿ����***%
I(:,:,1)=ccd1;
I(:,:,2)=ccd2;
I(:,:,3)=ccd3;
I(:,:,4)=ccd4;
I(:,:,5)=ccd5;

%***���������ֵ��ʼ��***%
S=zeros(m,n);
C=zeros(m,n);
sindelta=zeros(1,J);
cosdelta=zeros(1,J);
tandelta=zeros(1,J);
lamda=570;%��λnm
phi=zeros(m,n);%��������Ӧ����λ
% fs=phi*lamda/(2*pi);%��face shape����
kappa=1*10^(-3);
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

fs=phi*lamda/(2*pi);
mesh(phi);%�����άͼ