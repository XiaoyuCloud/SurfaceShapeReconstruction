% %�Ϻ����ݣ�19��41��42
clear all;
clc;

%***�趨һЩ��ʼֵ***%
J=9;
xun=[1,2,3,4,5,6,7,8,9];%ѡ���ļ�������ͼ���м���
m1=1;%m1��m2��m3��m4����ѡȡ����������ͼ���и���ͼ���ڵ�����
m2=200;
n1=1;
n2=200;
m=m2-m1+1;
n=n2-n1+1;
mn=110;
mcentre=round(m/2);
ncentre=round(n/2);
Mc=mcentre;
Nc=600;
I=zeros(m,n,J);%IΪ����ͼ����Ӧ�ĻҶ�ͼ�����ｫ���ʼ��
kappa=1*10^(-4);%��kappa�������о�
quyu1=81;%ѡ������������������XY��Ϊ��X�������Y�β������õ�������
quyu2=81;
Sumnum=quyu2-quyu1+1;
CharacterSum=zeros(Sumnum*5,9);%5*3
number=1;
lamda=570;%��λnm

refer=1;
meaed=1;
psimax=0.5*pi;
hphimax=20*pi;
hdelta=hphimax/n;
R=0.5*psimax+Nc^2/(2*psimax);
%***����һ������***%
Hphi=zeros(m,n);
phili=zeros(m,n);
deltali=[0,0.1,1*pi/5,2.5*pi/5,3*pi/5,4.2*pi/5,5.1*pi/5,6.5*pi/5,7*pi/5];
rd=500;
rdnum=3;
for q=1:m
    for p=1:n
        Hphi(p,q)=R-(R^2-((p-Mc)^2+(q+rd-Nc)^2))^(0.5);
    end;
    Hphi(:,q)=Hphi(:,q)+(q-1)*hdelta;
end;
% Hphi=Hphi+0.3*rand(m,n);
H=0.5*Hphi/(2*pi)*lamda;
deltaHphi=pi/5;
deltaH=0.5*deltaHphi/(2*pi)*lamda;%ÿ�βο���λ�ƶ���deltah��

%***��ȡ����ĸ���ͼ***%
for quyu=quyu1:quyu2
    jn=1;%JΪ����Ĵ�����
    ccdnum=zeros(m,n,J);
    cohnum=zeros(J,n);
    for j=1:9
        for p=1:m
            for q=1:n
%                 ccdnum(p,q,j)=abs(meaed+refer*cos(2*(H(p,q)-(j-1)*deltaH)*2*pi/lamda));%9������  
%                 phili(p,q)=2*H(p,q)*2*pi/lamda;
                ccdnum(p,q,j)=abs(meaed+refer*cos(Hphi(p,q)-deltali(j)));%9������  
                phili(p,q)=Hphi(p,q);
            end;
        end;
    end;
    for j=1:9
        I(:,:,jn)=ccdnum(:,:,jn);
        jn=jn+1;
    end;
    ccd1part=I(:,:,1);%����������벻���޸���
    figure(11)
    for j=1:9
        imshow(I(:,:,j));
    end;

    %***�趨�ο���λ�ƶ��ĳ�ʼֵ***%
    delta=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
%     delta=[0,0.1,1*pi/5,2.5*pi/5,3*pi/5,4.2*pi/5,5.1*pi/5,6.5*pi/5,7*pi/5];
      %�����ʼ����λ�ƶ���,��λ����
      %��ʼֵҪ�����˵�ӽ���ֵ���������ս��Ҳ���ܻ���ֽϴ�ƫ��
    deltachu=delta;
    
    %***��ʼ��***%
    deltak=zeros(1,J);%�����洢��һ��������õ���λ�ƶ���,��λ����
    deltaD=abs(delta-deltak);%����֮�����λ�仯��,��λ����
    psi=zeros(m,n,J);%psi�����ʼ��

    %***���������ֵ��ʼ��***%
    S=zeros(m,n);
    C=zeros(m,n);
    sindelta=zeros(1,J);
    cosdelta=zeros(1,J);
    tandelta=zeros(1,J);
    phi=zeros(m,n);%������������Ӧ����λ
    phichu=zeros(m,n);%atanֱ�ӵõ���phi
    phih=zeros(m,n);
    num=0;%��������

    %***�Ĳ�����������ݳ�ʼ��***%
    phifour=zeros(m,n);%�Ĳ����������õ���������Ӧ����λ
    sinfour=zeros(m,n);
    cosfour=zeros(m,n);

    %***����psi����***%
    for j=1:J
        psi(:,:,j)=I(:,:,j)-I(:,:,1);%psi����ֵ
    end;

    %***��С������ϵ����㷨***%
    while deltaD(1)>=kappa||deltaD(2)>=kappa||deltaD(3)>=kappa||deltaD(4)>=kappa||deltaD(5)>=kappa||deltaD(6)>=kappa||deltaD(7)>=kappa||deltaD(8)>=kappa||deltaD(9)>=kappa
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

    %***������Ӧ��phi***%
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
    %***��λ���***%
    for numsin=1:13 %�������λ���й�
        %mcentre��һ�к���
        p=mcentre;
        K=0;
        qsingul=zeros(1,10);
        cb=1;
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

        %����
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
    
    %***��ȥƫ��***%
    fs=0.5*phi*lamda/(2*pi);%��face shape
    deltabi=[deltali;deltachu;delta];

    %***�Ĳ������㷨***%
    deltafour=deltachu;%���ۣ�����ʼ��deltali������С������ϵõ���delta�����Ĳ����򷨼��������С���˷����������֮��,���ϵ����0.999���ϲűȽϺá�����������ϴ����ϵ���������½��ܶ࣬���Դ����ϵ�������ܿ�����ʲô
                      %deltafour��ʾ�Ĳ��������õĲο���λ
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

    %***��ȥƫ��***%
    phifourmin=min(phifour(:));
    phifour=phifour-phifourmin;

    picphi=figure(1);
    mesh(phi);%�����άͼ
%     axis([1 m 1 n -0.5 2*pi]);
    title('��С������ϵ����㷨');
    picphifour=figure(2);
    mesh(phifour);%�����άͼ
%     axis([1 m 1 n -0.5 2*pi]);
    title('�Ų������㷨');
    namephi=strcat('phi',num2str(quyu),'.jpg');
    namephifour=strcat('phifour',num2str(quyu),'.jpg');
    saveas(picphi,namephi);
    saveas(picphifour,namephifour);
    
    figure(100)
    for j=1:9
        subplot(2,5,j)
        imshow(I(:,:,j));
    end;
end;
nameCS=strcat(num2str(quyu1),'_',num2str(quyu2),'.xlsx');
xlswrite(nameCS,CharacterSum);


%����Ϊ��������
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
Z=phili(:);
phipolyli=fit([X,Y],Z,'poly11');
phipolyliZ=phipolyli(X,Y)';
phiplaneli=zeros(m,n);
phiplanelixie=zeros(m,n);
phiplanelixieli=zeros(m,n);
tiduquyu=(max(phili(:))-min(phili(:)))/(n-1);
for p=1:m
        phiplaneli(:,p)=(phili(:,p)-(phipolyliZ(1,((p-1)*m+1):p*m))');
        phiplanelixieli(:,p)=phili(:,p)-(p-1)*hdelta;
        phiplanelixie(:,p)=(phili(:,p)-(p-1)*tiduquyu);
end;

figure(101)
mesh(phiplaneli/(2*pi));
% axis([1 m 1 n -0.5 4*pi]);
title('��ȥ����γ�б��ƫ�ú�Ľ��');
figure(102)
mesh(phiplanelixie/(2*pi));
% axis([1 m 1 n -0.5 4*pi]);
title('��ȥֱ���γ�б��ƫ�ú�Ľ��');
figure(103)
mesh(phiplanelixieli/(2*pi));
% axis([1 m 1 n 0 2]);
title('���θ�ԭ���');
% title('��������');
figure(104)
plot(phipolyli,[X,Y],Z);
title('��ϵĽ��');

%����Ϊ�㷨�������������
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
phiplanexieli=zeros(m,n);
tiduquyu=(max(phili(:))-min(phili(:)))/(n-1);
for p=1:m
        phiplane(:,p)=(phi(:,p)-(phipolyZ(1,((p-1)*m+1):p*m))');
        phiplanexieli(:,p)=phi(:,p)-(p-1)*hdelta;%������-����ƫ��
        phiplanexie(:,p)=(phi(:,p)-(p-1)*tiduquyu);
end;
phiplanexielimin=min(phiplanexieli(:));
phiplanexieli=phiplanexieli-phiplanexielimin;

figure(3)
mesh(phiplane/(2*pi));
% axis([1 m 1 n -0.5 4*pi]);
title('��ȥ����γ�б��ƫ�ú�Ľ��');
figure(4)
mesh(phiplanexie/(2*pi));
% axis([1 m 1 n -0.5 4*pi]);
title('��ȥֱ���γ�б��ƫ�ú�Ľ��');
figure(5)
mesh(phiplanexieli/(2*pi));
% axis([1 m 1 n -0.5 4*pi]);
axis([1 m 1 n 0 0.5]);
title('���θ�ԭ���');
% title('��ȥ����б��ƫ�ú�Ľ��');
figure(6)
plot(phipoly,[X,Y],Z);
title('��ϵĽ��');

phiplane=phiplanexieli;
Rdelta=corrcoef(delta,deltafour);
Rphi=corrcoef(phi,phifour);%�����㷨�������㷨������õ����ε����ϵ��
MEANphiplane=mean(phiplane(:));
MEANphifour=mean(phifour(:));
RMSphiplane=rms(phiplane(:));
RMSphifour=rms(phifour(:));
MAXMINphiplane=max(phiplane(:))-min(phiplane(:));
MAXMINphifour=max(phifour(:))-min(phifour(:));
VARphiplane=var(phiplane(:));
VARphifour=var(phifour(:));
character=[MEANphiplane,MAXMINphiplane,RMSphiplane,VARphiplane,0,0,0,0,0;MEANphifour,MAXMINphifour,RMSphifour,VARphifour,0,0,0,0,0];
CharacterSum(number:number+4,:)=[deltabi;character];
number=number+6;

namematphiplanelixieli=strcat('lixieli',num2str(rdnum),'.mat');
save(namematphiplanelixieli,'phiplanelixieli');
