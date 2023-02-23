% x =1:m;  %% x,y取值范围
% y =1:n;
% [xX, yY] = meshgrid(x,y);  %% 生成一个二维的取值范围
% [mM, nN] = size(xX);
% 
% x = reshape(xX, mM*nN, 1);  %% 把矩阵转化为向量
% y = reshape(yY, mM*nN, 1);

% 待拟合方程：F = z^2 = (-c^2/a^2*x^2) + (c^2/a^2*2*x1*x) + (- c^2/b^2*y^2) +
%                      (c^2/b^2*2*y1*y) + (2*z1*z) +
%                      (-c^2/a^2*x1^2 - c^2/b^2*y1^2 - z1^2 + c^2)
% x,y,z 均要先转化为列向量
% k(1) = -c^2/a^2  由k值就可求出椭圆所有参数！！！
% k(2) = c^2/a^2*2*x1
% k(3) = - c^2/b^2
% k(4) = c^2/b^2*2*y1
% k(5) = 2*z1
% k(6) = -c^2/a^2*x1^2 - c^2/b^2*y1^2 - z1^2 + c^2

z=phiquyuplane(:);
xdata = [X,Y,z];  %% 将 x，y 数据按列组合到 xdata
ydata = z.^2;  %% 先把 z 值平方，再进行拟合

k0 = ones(1,6);  %% k 的运行初值，不会影响最终结果

F = @(k,xdata) k(1)*xdata(:,1).^2 + k(2)*xdata(:,1) + k(3)*xdata(:,2).^2 + k(4)*xdata(:,2) + k(5)*xdata(:,3) + k(6);
[k,resnorm]=lsqcurvefit(F,k0,xdata,ydata);

% step2：椭圆参数求解
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
disp('a轴:');
disp(a);
disp('b轴:');
disp(b);
disp('c轴:');
disp(c);

phiquyutuo=zeros(m,n);
for p=1:m
    for q=1:n
        phiquyutuo(p,q)=sqrt((1-(p-x1)^2/a^2 -(q-y1)^2/b^2)*c^2)+z1;
    end;
end;
figure(14)
mesh(phiquyutuo);