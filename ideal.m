%�Ϻ����ݣ�19��41��42
clc;

%***�趨һЩ��ʼֵ***%
J=9;
xun=[1,2,3,4,5,6,7,8,9];%ѡ���ļ�������ͼ���м���
m1=494;%m1��m2��m3��m4����ѡȡ����������ͼ���и���ͼ���ڵ�����
m2=648;
n1=468;
n2=612;
m=m2-m1+1;
n=n2-n1+1;
I=zeros(m,n,J);%IΪ����ͼ����Ӧ�ĻҶ�ͼ�����ｫ���ʼ��
kappa=1*10^(-3);%��kappa�������о�
quyu1=81;%ѡ������������������XY��Ϊ��X�������Y�β������õ�������
quyu2=81;
Sumnum=quyu2-quyu1+1;
CharacterSum=zeros(Sumnum*6,9);%5*3
number=1;

ff=n;
r=(0:ff-1)/ff;
A=1;
f=1;

%***��ȡ����ĸ���ͼ***%
for quyu=quyu1:quyu2
    jn=1;%JΪ����Ĵ�����
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
            ccdnum(mm,:,j)=cohnum(j,:);%ccd����ͼ
        end;
    end;
    for j=1:9
        I(:,:,jn)=ccdnum(:,:,jn);
        jn=jn+1;
    end;
%     imshow(I(:,:,3));

    %***�趨�ο���λ�ƶ��ĳ�ʼֵ***%
%     delta=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
    delta=[0,0.1,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5];
      %�����ʼ����λ�ƶ���,��λ����
      %��ʼֵҪ�����˵�ӽ���ֵ���������ս��Ҳ���ܻ���ֽϴ�ƫ��
    deltachu=delta;
    deltali=[0,1*pi/5,2*pi/5,3*pi/5,4*pi/5,5*pi/5,6*pi/5,7*pi/5,8*pi/5];
    
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
    lamda=570;%��λnm
    phi=zeros(m,n);%������������Ӧ����λ
    phichu=zeros(m,n);%atanֱ�ӵõ���phi
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
                phi(p,q)=atan(S(p,q)/C(p,q))+2*pi;
            elseif (S(p,q)>0&&C(p,q)==0)
                phi(p,q)=pi/2;
            elseif (S(p,q)<0&&C(p,q)==0)
                phi(p,q)=3*pi/2;
            end;
        end;
    end;

    %***��ȥƫ��***%
%     phimin=min(phi(:));
%     phi=phi-phimin;
    fs=0.5*phi*lamda/(2*pi);%��face shape
    deltabi=[deltali;deltachu;delta;deltali-delta];

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
                    phifour(p,q)=atan(sinfour(p,q)/cosfour(p,q))+2*pi;
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
    title('��С������ϵ����㷨');
%     axis([1 m 1 n 0 2*pi]);
    picphifour=figure(2);
    mesh(phifour);%�����άͼ
    title('�Ų������㷨');
    namephi=strcat('phi',num2str(quyu),'.jpg');
    namephifour=strcat('phifour',num2str(quyu),'.jpg');
    saveas(picphi,namephi);
    saveas(picphifour,namephifour);
    
    Rdelta=corrcoef(delta,deltafour);
    Rphi=corrcoef(phi,phifour);%�����㷨�������㷨������õ����ε����ϵ��
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
