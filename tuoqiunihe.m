% x =1:m;  %% x,yȡֵ��Χ
% y =1:n;
% [xX, yY] = meshgrid(x,y);  %% ����һ����ά��ȡֵ��Χ
% [mM, nN] = size(xX);
% 
% x = reshape(xX, mM*nN, 1);  %% �Ѿ���ת��Ϊ����
% y = reshape(yY, mM*nN, 1);

% ����Ϸ��̣�F = z^2 = (-c^2/a^2*x^2) + (c^2/a^2*2*x1*x) + (- c^2/b^2*y^2) +
%                      (c^2/b^2*2*y1*y) + (2*z1*z) +
%                      (-c^2/a^2*x1^2 - c^2/b^2*y1^2 - z1^2 + c^2)
% x,y,z ��Ҫ��ת��Ϊ������
% k(1) = -c^2/a^2  ��kֵ�Ϳ������Բ���в���������
% k(2) = c^2/a^2*2*x1
% k(3) = - c^2/b^2
% k(4) = c^2/b^2*2*y1
% k(5) = 2*z1
% k(6) = -c^2/a^2*x1^2 - c^2/b^2*y1^2 - z1^2 + c^2

z=phiquyuplane(:);
xdata = [X,Y,z];  %% �� x��y ���ݰ�����ϵ� xdata
ydata = z.^2;  %% �Ȱ� z ֵƽ�����ٽ������

k0 = ones(1,6);  %% k �����г�ֵ������Ӱ�����ս��

F = @(k,xdata) k(1)*xdata(:,1).^2 + k(2)*xdata(:,1) + k(3)*xdata(:,2).^2 + k(4)*xdata(:,2) + k(5)*xdata(:,3) + k(6);
[k,resnorm]=lsqcurvefit(F,k0,xdata,ydata);

% step2����Բ�������
x1 = -k(2)/k(1)/2;
y1 = -k(4)/k(3)/2;
z1 = k(5)/2;
c = abs(sqrt((z1^2 + k(6))/(1 - 1/a^2*x1^2 - 1/b^2*y1^2)));
a = c/abs(sqrt(-k(1)));
b = c/abs(sqrt(-k(3)));

disp('x1:');
disp(x1);
disp('y1:');
disp(y1);
disp('z1:');
disp(z1);
disp('a��:');
disp(a);
disp('b��:');
disp(b);
disp('c��:');
disp(c);

phiquyutuo=zeros(m,n);
for p=1:m
    for q=1:n
        phiquyutuo(p,q)=sqrt((1-(p-x1)^2/a^2 -(q-y1)^2/b^2)*c^2)+z1;
    end;
end;
figure(14)
mesh(phiquyutuo);