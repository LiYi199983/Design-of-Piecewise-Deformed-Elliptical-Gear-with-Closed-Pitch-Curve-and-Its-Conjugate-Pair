clear;clc;close all
k1=0.1;%主动轮偏心率
XN1=1;      %主动轮阶数（不可变动）
XN2=1;      %从动轮阶数
mn=4;%法向模数
Bc=13.8717*pi/180;%螺旋角
m11=1.6       ;%可不为整数
m12=2;            %0.5;可不为整数
m13=1/(3-1/m11-1/m12);%m11*m12/(3*m11*m12-m11-m12);            %0.5;可不为整数
z1=44;%主动轮齿数【推荐值】%%%%%%%%高阶椭圆齿轮齿数为其阶数和某一奇数的乘积【mod(z1*XN2,XN1)必须等于0】
w1= 0.6283;%主动轮角速度

i=1;  
zz1=z1;
while mod(z1*XN2,XN1)>0 
z1=zz1+round((i*(-1)^i)/2);%双向寻找被XN1整除的z1*XN2（z1）
i=i+1;
end
z2=z1*XN2/XN1;  %从动轮齿数%从动轮齿数%从动轮齿数
disp(['z1(主动轮齿数）=' num2str(z1)]);
disp(['z2(从动轮齿数）=' num2str(z2)]);
mt=mn/cos(Bc);%端面模数mm
an=20*pi/180;%刀具齿形角
at=atan(tan(an)/cos(Bc));%端面压力角

DN=3 ;      %分段数(不可变动)
L1=pi*mt*z1;%整圆周长
L2=L1*XN2/XN1;%整圆周长

tttt=1.5;%%%%%%%%%%%%%%%%%%%%%%%%%%保证1.5秒加工完成一圈%%%%%%是“jiange”倍数
jiange=0.0015; %0.003;%%%%%%%%%%%%%%%%%%%%%%%%%间隔时间%%%%%×3
Temjiange=0.0001;
han=1;
cn=0.25;
hat=han*cos(Bc);
ct=cn*cos(Bc);
ha1=han*mn;
ha2=han*mn;
hf1=(han+cn)*mn;
hf2=(han+cn)*mn;


%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S
 syms sitaF1
 syms sitaF2
 syms sitaF3
 rr11=1/(1-k1*cos(XN1*m11*sitaF1));%极径
 drr11dsita=diff(rr11, sitaF1,1);
 fun4=inline(sqrt(rr11^2+drr11dsita.^2)); 
 sss1=quadl(fun4,0,2*pi/(DN*XN1*m11));
   
 rr22=1/(1-k1*cos(XN1*m12*(sitaF2-2*pi/(DN*XN1*m11))+1*2*pi/DN));%极径
 drr22dsita=diff(rr22, sitaF2,1);
 fun5=inline(sqrt(rr22^2+drr22dsita.^2)); 
 sss2=quadl(fun5,2*pi/(DN*XN1*m11),2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11));
  
 rr33=1/(1-k1*cos(XN1*m13*(sitaF3-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN));%极径
 drr33dsita=diff(rr33, sitaF3,1);
 fun6=inline(sqrt(rr33^2+drr33dsita.^2)); 
 sss3=quadl(fun6,2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11),2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11)+2*pi/(DN*XN1*m13));
 sss=XN1*(sss1+sss2+sss3);
%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S%计算弧长S
A1=L1/(sss*(1-k1^2));% %数值积分法计算主动轮半长轴
disp(['A1(主动轮半长轴）=' num2str(A1)]);

%根据封闭条件求中心距
A1A2=zhxinjunJ3D(k1,XN1,XN2,A1,DN,m11,m12,m13);% %黄金分割法计算半长轴%黄金分割法计算半长轴
disp(['A1A2(中心距）=' num2str(A1A2)]);
A2=A1A2-A1;
 
p1=A1*(1-k1^2);
axis equal 
hold on

for i=1:1:tttt/jiange;
 sita1(i)=i*((2*pi)/tttt)*jiange;     
   
for jjj=1:1:XN1    %jjj周期号
if (sita1(i)-2*pi*(jjj-1)/XN1)>0 &&(sita1(i)-2*pi*(jjj-1)/XN1)<=2*pi/(DN*XN1*m11);  
  r1(i)=p1/(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)));%极径
  i12(i)=(A1A2-p1-A1A2*k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))/p1;%传动比%非圆齿轮
  u1(i)=acos(-(k1*XN1*m11*sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))/((k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))- 1)^2+k1^2*(XN1*m11)^2*(sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2 )^(1/2) );
  a12(i)=(u1(i)+at-pi/2)*180/pi;%主动轮的压力角
  
  qlbj1(i)=p1*((1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2+k1^2*(XN1*m11)^2*(sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2)^1.5/( (1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^3*(1+ k1*((XN1*m11)^2-1)* cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)) ));%主动轮的曲率半径11
  minqlbj1(i)=qlbj1(i);   
  Fgun1(i)= 1+((XN1*m11)^2-1)*k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1));%%主动轮节曲线上无内凹部分的条件
  if ((XN1*m11)^2-1)*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))<0
  KK11max(i)=1/((1-(XN1*m11)^2)*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)));%主动轮非内凹最大偏心率
  else 
  KK11max(i)=1;%此时应该是无穷大，意即任意k1均满足。  
  end
  FBZH1=i;
  
  dr11dsita1(i)=-p1*k1*(XN1*m11)* sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))/(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2;%主动轮一阶导数
  d2r11dsita12(i)=A1*k1*(XN1*m11)^2*(1-k1^2)*(2*k1*(sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2-cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))*(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))))/(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^3;%主动轮二阶导数
  dr21dsita2(i)=-dr11dsita1(i)*i12(i);%从动轮一阶导数
  d2r21dsita22(i)=-d2r11dsita12(i)*i12(i)^2-dr11dsita1(i)*(XN1*m11)*A1A2*k1*sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))*i12(i)/p1; %从动轮二阶导数  
% dr21dsita2o(i)=(XN1*m11)*k1*sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))*(A1A2-p1-A1A2*k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))/(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2;%从动轮一阶导数（另一种形式）
% d2r21dsita2o(i)=(XN1*m11)^2*k1*(A1A2-p1-A1A2*k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2*((cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))-k1*(sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2-k1)/(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))+k1*A1A2*(sin(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2/(A1A2-p1-A1A2*k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1))))/(p1*(1-k1*cos(XN1*m11*(sita1(i)-2*pi*(jjj-1)/XN1)))^2);%从动轮二阶导数（另一种形式） 
%  QQQ1(i)=dr21dsita2o(i)-dr21dsita2(i); 
%  QQQ2(i)=d2r21dsita2o(i)-d2r21dsita22(i);  
  
  
  
  
elseif sita1(i)-2*pi*(jjj-1)/XN1>2*pi/(DN*XN1*m11)&&(sita1(i)-2*pi*(jjj-1)/XN1)<=2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11);  
  r1(i)=p1/(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN));%极径
  i12(i)=(A1A2-p1-A1A2*k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))/p1;%传动比%非圆齿轮
  u1(i)=acos(-(k1*XN1*m12*sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))/((k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)- 1)^2+k1^2*(XN1*m12)^2*(sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2 )^(1/2) );
  a12(i)=(u1(i)+at-pi/2)*180/pi;%主动轮的压力角
  
  qlbj1(i)=p1*((1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2+k1^2*(XN1*m12)^2*(sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2)^1.5/( (1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^3*(1+ k1*((XN1*m12)^2-1)* cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN) ));%主动轮的曲率半径12
  minqlbj2(i)=qlbj1(i);   
  Fgun1(i)= 1+((XN1*m12)^2-1)*k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN);%%主动轮节曲线上无内凹部分的条件
  if ((XN1*m12)^2-1)*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)<0
  KK12max(i)=1/((1-(XN1*m12)^2)*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN));%主动轮非内凹最大偏心率
  else 
  KK12max(i)=1;%此时应该是无穷大，意即任意k1均满足。  
  end
  FBZH2=i;
  
dr11dsita1(i)=-p1*k1*(XN1*m12)* sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)/ (1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2;%主动轮一阶导数
d2r11dsita12(i)=A1*k1*(XN1*m12)^2*(1-k1^2)*(2*k1*(sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2-cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)*(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)))/(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^3;%主动轮二阶导数
dr21dsita2(i)=-dr11dsita1(i)*i12(i);%从动轮一阶导数
d2r21dsita22(i)=-d2r11dsita12(i)*i12(i)^2-dr11dsita1(i)*(XN1*m12)*A1A2*k1*sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)*i12(i)/p1;%从动轮二阶导数     
% dr21dsita2o(i)=(XN1*m12)*k1*sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)*(A1A2-p1-A1A2*k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))/(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2;%从动轮一阶导数（另一种形式）
% d2r21dsita2o(i)=(XN1*m12)^2*k1*(A1A2-p1-A1A2*k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2*((cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)-k1*(sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2-k1)/(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))+k1*A1A2*(sin(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2/(A1A2-p1-A1A2*k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN)))/(p1*(1-k1*cos(XN1*m12*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11))+2*pi/DN))^2);%从动轮二阶导数（另一种形式）  
%  QQQ1(i)=dr21dsita2o(i)-dr21dsita2(i); 
%  QQQ2(i)=d2r21dsita2o(i)-d2r21dsita22(i);  

elseif sita1(i)-2*pi*(jjj-1)/XN1>2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11)&&(sita1(i)-2*pi*(jjj-1)/XN1)<=2*pi/(DN*XN1*m13)+2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m11);  
  r1(i)=p1/(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN));%极径
  i12(i)=(A1A2-p1-A1A2*k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))/p1;%传动比%非圆齿轮
  u1(i)=acos(-(k1*XN1*m13*sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))/((k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)- 1)^2+k1^2*(XN1*m13)^2*(sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2 )^(1/2) );
  a12(i)=(u1(i)+at-pi/2)*180/pi;%主动轮的压力角
  
  qlbj1(i)=p1*((1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2+k1^2*(XN1*m13)^2*(sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2)^1.5/( (1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^3*(1+ k1*((XN1*m13)^2-1)* cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN) ));%主动轮的曲率半径11
  minqlbj3(i)=qlbj1(i);   
  Fgun1(i)= 1+((XN1*m13)^2-1)*k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN);%%主动轮节曲线上无内凹部分的条件
  if ((XN1*m13)^2-1)*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)<0
  KK13max(i)=1/((1-(XN1*m13)^2)*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN));%主动轮非内凹最大偏心率
   else 
  KK13max(i)=1;%此时应该是无穷大，意即任意k1均满足。  
  end
  FBZH3=i;
  
dr11dsita1(i)=-p1*k1*(XN1*m13)* sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)/ (1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2;%主动轮一阶导数
d2r11dsita12(i)=A1*k1*(XN1*m13)^2*(1-k1^2)*(2*k1*(sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2-cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)*(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)))/(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^3;%主动轮二阶导数
dr21dsita2(i)=-dr11dsita1(i)*i12(i);%从动轮一阶导数
d2r21dsita22(i)=-d2r11dsita12(i)*i12(i)^2-dr11dsita1(i)*(XN1*m13)*A1A2*k1*sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)*i12(i)/p1;%从动轮二阶导数     
%  dr21dsita2o(i)=(XN1*m13)*k1*sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)*(A1A2-p1-A1A2*k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))/(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2;%从动轮一阶导数（另一种形式）
%  d2r21dsita2o(i)=(XN1*m13)^2*k1*(A1A2-p1-A1A2*k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2*((cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)-k1*(sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2-k1)/(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))+k1*A1A2*(sin(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2/(A1A2-p1-A1A2*k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)))/(p1*(1-k1*cos(XN1*m13*(sita1(i)-2*pi*(jjj-1)/XN1-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN))^2);%从动轮二阶导数（另一种形式）  
%   QQQ1(i)=dr21dsita2o(i)-dr21dsita2(i); 
%   QQQ2(i)=d2r21dsita2o(i)-d2r21dsita22(i);  
end  
end

 if sita1(i)>=0 && sita1(i)<=2*pi/(DN*XN1*m11);  
   for j=1:1:XN2 
 TemF1=fix(tttt*(j-1)/(XN1*jiange));   
 Temsita1=0:Temjiange:sita1(i);
 TemJFX1=p1./(A1A2-p1-A1A2*k1*cos(XN1*m11*Temsita1));  
 JCHsita11=trapz(TemJFX1)*Temjiange;%本循环是单次值；后面用的是最终值
 sita2(TemF1+i)=2*pi*(j-1)/(XN2)+JCHsita11;   
 r2(TemF1+i)= A1A2-r1(i);
  qlbj2(TemF1+i)=(r2(TemF1+i)^2+(dr21dsita2(i))^2)^1.5/(r2(TemF1+i)^2+2*(dr21dsita2(i))^2-r2(TemF1+i)*d2r21dsita22(i));
  w2(TemF1+i)=w1/i12(i); %从动轮角速度
   di12dsita1(i)=A1A2*k1*XN1*m11*sin(XN1*m11*(sita1(i)))/p1;
   alfa2(TemF1+i)=-w1^2* di12dsita1(i)/ i12(i)^2;%从动轮角加速度

   end   
 elseif sita1(i)>=2*pi/(DN*XN1*m11) && sita1(i)<=2*pi/(DN*XN1*m11)+2*pi/(DN*XN1*m12);  
 for j=1:1:XN2 
 TemF1=fix(tttt*(j-1)/(XN1*jiange));  
 Temsita2=2*pi/(XN1*DN*m11):Temjiange:sita1(i);
 TemJFX2=p1./(A1A2-p1-A1A2*k1*cos(XN1*m12*(Temsita2-2*pi/(XN1*DN*m11))+2*pi/DN));  
 JCHsita22=trapz(TemJFX2)*Temjiange;%本循环是单次值；后面用的是最终值
 sita2(TemF1+i)=2*pi*(j-1)/(XN2)+JCHsita11+JCHsita22;   
 r2(TemF1+i)= A1A2-r1(i);
 qlbj2(TemF1+i)=(r2(TemF1+i)^2+(dr21dsita2(i))^2)^1.5/(r2(TemF1+i)^2+2*(dr21dsita2(i))^2-r2(TemF1+i)*d2r21dsita22(i));
  w2(TemF1+i)=w1/i12(i); %从动轮角速度
    di12dsita1(i)=A1A2*k1*XN1*m12*sin(XN1*m12*(sita1(i)-2*pi/(DN*XN1*m11))+2*pi/DN)/p1;
   alfa2(TemF1+i)=-w1^2* di12dsita1(i)/ i12(i)^2;%从动轮角加速度
 end 
 
  elseif sita1(i)>=2*pi/(DN*XN1*m11)+2*pi/(DN*XN1*m12) && sita1(i)<=2*pi/(DN*XN1*m11)+2*pi/(DN*XN1*m12)+2*pi/(DN*XN1*m13);  
 for j=1:1:XN2 
 TemF1=fix(tttt*(j-1)/(XN1*jiange));  
 Temsita3=2*pi/(XN1*DN*m12)+2*pi/(XN1*DN*m11):Temjiange:sita1(i);
 TemJFX3=p1./(A1A2-p1-A1A2*k1*cos(XN1*m13*(Temsita3-2*pi/(XN1*DN*m11)-2*pi/(XN1*DN*m12))+2*2*pi/DN));   
 sita2(TemF1+i)=2*pi*(j-1)/(XN2)+trapz(TemJFX3)*Temjiange+JCHsita11+JCHsita22;   
 r2(TemF1+i)= A1A2-r1(i);  
 qlbj2(TemF1+i)=(r2(TemF1+i)^2+(dr21dsita2(i))^2)^1.5/(r2(TemF1+i)^2+2*(dr21dsita2(i))^2-r2(TemF1+i)*d2r21dsita22(i));
 w2(TemF1+i)=w1/i12(i); %从动轮角速度
   di12dsita1(i)=A1A2*k1* XN1*m13*sin(XN1*m13*(sita1(i)-2*pi/(DN*XN1*m11)-2*pi/(DN*XN1*m12))+2*2*pi/DN)/p1;
   alfa2(TemF1+i)=-w1^2* di12dsita1(i)/ i12(i)^2;%从动轮角加速度   
 end 
 end   
 
  for jjj=1:1:XN1
  TemF1=fix(tttt*(jjj-1)/(XN1*jiange));    
  if sita1(i)>=0 && sita1(i)<=2*pi/XN1;  
  UU1(TemF1+i)=sqrt((abs(qlbj1(i))+ha1)^2- abs(qlbj1(i))^2*cos(Bc)^2/(tan(an)^2+cos(Bc)^2))- abs(qlbj1(i))* sqrt(tan(an)^2/(tan(an)^2+cos(Bc)^2));
  UU2(TemF1+i)=sqrt((abs(qlbj2(i))+ha2)^2- abs(qlbj2(i))^2*cos(Bc)^2/(tan(an)^2+cos(Bc)^2))- abs(qlbj2(i))* sqrt(tan(an)^2/(tan(an)^2+cos(Bc)^2));
  Chhedu(TemF1+i)=( UU1(TemF1+i)+ UU2(TemF1+i))* sqrt(tan(an)^2+cos(Bc)^2)/(pi*mn);%啮合重合度%啮合重合度%啮合重合度%啮合重合度%啮合重合度
  sita1sita1(TemF1+i)=(jjj-1)*2*pi/XN1+sita1(i);
  end
  end

  
 
%  u1(i)=pi+u1(i)
Ra1(i)=sqrt(ha1^2+r1(i)^2+2* r1(i)* ha1*sin( u1(i)));%主动轮齿顶圆极径
if u1(i)>=pi/2
Sitaa1(i)=sita1(i)+asin(ha1*sin(u1(i))/ Ra1(i));%主动轮齿顶圆极角
else
Sitaa1(i)=sita1(i)-asin(ha1*sin(u1(i))/ Ra1(i));%主动轮齿顶圆极角
end

Rf1(i)=sqrt(hf1^2+r1(i)^2-2* r1(i)* hf1*sin( u1(i)));%主动轮齿根圆极径
if u1(i)>=pi/2
Sitaf1(i)=sita1(i)-asin(hf1*sin(u1(i))/ Rf1(i));%主动轮齿顶圆极角
else
Sitaf1(i)=sita1(i)+asin(hf1*sin(u1(i))/ Rf1(i));%主动轮齿顶圆极角
end
    
end

%齿顶曲线与齿根曲线自交问题及解决%齿顶曲线与齿根曲线自交问题及解决%齿顶曲线与齿根曲线自交问题及解决
FJDRa=0; 
FJDRf=0; 
for i=1:1:tttt/jiange;
Fi1=i-1;
if Fi1<=0
  Fi1=Fi1+fix(tttt/jiange);  
end    
if Sitaa1(i)<Sitaa1(Fi1)
  FJDRa=1;
   j=0;
end
if  FJDRa==1
     j=j+1;
Fi2=i-1-2*j;
if Fi2<=0
  Fi2=Fi2+fix(tttt/jiange);  
end        
  if Sitaa1(Fi2)<=Sitaa1(i)   
   for k=0:1:1+2*j 
Fi3=i-k;
if Fi3<=0
  Fi3=Fi3+fix(tttt/jiange);  
end  
    Sitaa1(Fi3)=Sitaa1(i);
    Ra1(Fi3)=Ra1(i);
   end
   FJDRa=0; 
end     
end  
end

for i=1:1:tttt/jiange;
Fi1=i-1;
if Fi1<=0
  Fi1=Fi1+fix(tttt/jiange);  
%  Sitaf1(Fi1)=Sitaf1(Fi1)-2*pi;
end   
if Sitaf1(i)<Sitaf1(Fi1)
  FJDRf=1;
   jj=0;
end
if  FJDRf==1
     jj=jj+1; 
   Fi2=i-1-2*jj;
if Fi2<=0
  Fi2=Fi2+fix(tttt/jiange);  
  Sitaf1(Fi2)=Sitaf1(Fi2)-2*pi;
end   
  if Sitaf1(Fi2)<=Sitaf1(i)   
   for k=0:1:1+2*jj 
  Fi3=i-k;
if Fi3<=0
  Fi3=Fi3+fix(tttt/jiange); 
end       
    Sitaf1(Fi3)=Sitaf1(i);
    Rf1(Fi3)=Rf1(i);
   end
   FJDRf=0; 
end     
end    
end    
%齿顶曲线与齿根曲线自交问题及解决%齿顶曲线与齿根曲线自交问题及解决%齿顶曲线与齿根曲线自交问题及解决

k11max=min(KK11max(FBZH3*(XN1-1)/XN1+1:FBZH1));%主动轮非内凹最大偏心率
k12max=min(KK12max(FBZH1+1:FBZH2));%主动轮非内凹最大偏心率
k13max=min(KK13max(FBZH2+1:FBZH3));%主动轮非内凹最大偏心率
disp(['k11max(主动轮非内凹最大偏心率）=' num2str(k11max) '主动轮实际偏心率k11(应小于前者）=' num2str(k1)]);
disp(['k12max(主动轮非内凹最大偏心率）=' num2str(k12max) '主动轮实际偏心率k12(应小于前者）=' num2str(k1)]);
disp(['k13max(主动轮非内凹最大偏心率）=' num2str(k13max) '主动轮实际偏心率k13(应小于前者）=' num2str(k1)]);

minqlbj1=minqlbj1(FBZH3*(XN1-1)/XN1+1:FBZH1);   
minqlbj2=minqlbj2(FBZH1+1:FBZH2);   
minqlbj3=minqlbj3(FBZH2+1:FBZH3);   
qlbj11min=min(minqlbj1);%主动轮最小曲率半径
qlbj12min=min(minqlbj2);%主动轮最小曲率半径
qlbj13min=min(minqlbj3);%主动轮最小曲率半径
disp(['qlbj11min(主动轮最小曲率半径）=' num2str(qlbj11min)]);
disp(['qlbj12min(主动轮最小曲率半径）=' num2str(qlbj12min)]);
disp(['qlbj13min(主动轮最小曲率半径）=' num2str(qlbj13min)]);

for i=1:1:tttt/jiange;
if qlbj1(i)<0
  qlbj1jdzhi(i)=100000;  
else
  qlbj1jdzhi(i)=qlbj1(i);    
end    
end
qlbj1min=min(qlbj1jdzhi);%主动轮最小曲率半径(实际有用的)



for i=1:1:tttt/jiange;
if qlbj2(i)<0
  qlbj2jdzhi(i)=100000;  
else
  qlbj2jdzhi(i)=qlbj2(i);    
end    
end
qlbj2min=min(qlbj2jdzhi);%从动轮最小曲率半径
disp(['qlbj2min(从动轮最小曲率半径）=' num2str(qlbj2min)]);

 
 mnmax1=qlbj1min*(tan(an))^2/(((cos(Bc))^2+(tan(an))^2)*han);%滚切主动轮不根切法面最大模数
 mnmax2=qlbj2min*(tan(an))^2/(((cos(Bc))^2+(tan(an))^2)*han);%滚切从动轮不根切法面最大模数
 disp(['mnmax1(滚切主动轮不根切法面最大模数）=' num2str(mnmax1) 'mn(实际法面模数）=' num2str(mn)]);
 disp(['mnmax2(滚切从动轮不根切法面最大模数）=' num2str(mnmax2) 'mn(实际法面模数）=' num2str(mn)]);
% 
z0min1=(mn^2*han^2*(cos(Bc))^3+ mn^2*han^2*cos(Bc)* (tan(an))^2- qlbj1min^2*(tan(an))^2*cos(Bc))/( qlbj1min* (tan(an))^2*mn- mn^2*han* (cos(Bc))^2- mn^2*han*(tan(an))^2);
z0min2=(mn^2*han^2*(cos(Bc))^3+ mn^2*han^2*cos(Bc)* (tan(an))^2- qlbj2min^2*(tan(an))^2*cos(Bc))/( qlbj2min* (tan(an))^2*mn- mn^2*han* (cos(Bc))^2- mn^2*han*(tan(an))^2);
 disp(['z0min1(插削主动轮不根切插齿刀最小齿数）=' num2str(z0min1)]);
 disp(['z0min2(插削从动轮不根切插齿刀最小齿数）=' num2str(z0min2)]);
% syms xxx
% solve(mt^2*hat*(xxx+hat)-min(qlbj11min,qlbj12min)*sin(at)^2*(min(qlbj11min,qlbj12min)+mt*xxx),'xxx') %数值解 
% solve(mt^2*hat*(xxx+hat)-min(qlbj21min,qlbj22min)*sin(at)^2*(min(qlbj21min,qlbj22min)+mt*xxx) , ' xxx ' ) %数值解


%polar(sita,r,'--r');%画主动轮节曲线
figure(1)
axis equal 
hold on
r1x=r1.*cos(sita1);
r1y=r1.*sin(sita1);
plot(r1x,r1y,'Color','b','LineStyle','-','linewidth',2);%画主动轮节曲线
% qlbj2=qlbj22;
r2x=r2.*cos(-sita2+pi)+A1A2;
r2y=r2.*sin(-sita2+pi);
ylabel('\ity\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 14);
xlabel('\itx\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 14);
set(gca,'XTick',-80:40:240,'Fontname', 'Times New Roman', 'Fontsize', 14)
set(gca,'YTick',-120:40:120,'Fontname', 'Times New Roman', 'Fontsize', 14)
plot(r2x,r2y,'Color','r','LineStyle','-.','linewidth',2);%画从动轮节曲线
plot(A1A2,0, 'b+');%画从动轮中心点
plot(0,0, 'r+');%画主动轮中心点
legend({'Driving Gear','Driven Gear'},'Fontsize', 14,'fontname','times new roman');

% title('齿轮副节曲线') 

figure(2)%传动比&主动轮的压力角
hold on
[AX,H1,H2] = plotyy(sita1,i12, sita1,a12,'plot');
set(AX(:),'Ycolor','k','Fontsize', 14); %设定两个Y轴的刻度字号
set(AX(1),'ylim',[0.85,1.35],'ytick',0.85:0.1:1.35,'FontName','Times New Roman');  %左Y轴的范围
set(AX(2),'ylim',[10,35],'ytick',10:5:35,'FontName','Times New Roman');  %右Y轴的范围
set(AX(1),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %左X轴的范围
set(AX(2),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %右X轴的范围
set(AX,'XTickLabel',{'0','0.5π','π','1.5π','2π'},'Fontsize', 14,'fontname','times new roman');%X轴标签
set(H1,'LineStyle','-','LineWidth',2,'color','b');%设置左坐标系曲线线型、标记及颜色
set(H2,'LineStyle','-.','LineWidth',2,'color','r');%设置右坐标系曲线线型、标记及颜色
%以下设置坐标系参数
set(AX(1),'XColor','k','YColor','b');%设置左坐标系颜色
set(AX(2),'XColor','k','YColor','r');%设置右坐标系颜色
set(get(AX(1),'Ylabel'),'string','Gear Ratio\it i','color','b','Fontsize', 14,'fontname','times new roman');%%左侧坐标轴label黑色字号字体
set(get(AX(2),'Ylabel'),'string','Pressure Angle\it \gamma\rm (^o )','color','r','Fontsize', 14,'fontname','times new roman');%%右侧坐标轴label字体黑色字号字体
xlabel('Polar Angle of Driving Gear {\itθ} (rad)','color','k','Fontsize', 14,'fontname','times new roman');%设置X轴标签
H1H=legend([H1,H2],{'Gear Ratio','Pressure Angle'},'Fontsize', 14,'FontName','Times New Roman');%HHHH1,'Δαc',

%legend('boxoff');  %legend去掉框
%box off;
hold off



figure(3)%主动轮的曲率半径&主动轮的节曲线上无内凹部分的条件
hold on
[AAX,HH1,HH2] = plotyy(sita1,qlbj1, sita1,Fgun1,'plot');
set(AAX(:),'Ycolor','k','Fontsize', 14); %设定两个Y轴的刻度字号
set(AAX(1),'ylim',[0,140],'ytick',0:20:140,'FontName','Times New Roman');  %左Y轴的范围
set(AAX(2),'ylim',[0.6,2.0],'ytick',0.6:0.2:2.0,'FontName','Times New Roman');  %右Y轴的范围
set(AAX(1),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %左X轴的范围
set(AAX(2),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %右X轴的范围
set(AAX,'XTickLabel',{'0','0.5π','π','1.5π','2π'},'Fontsize', 14,'fontname','times new roman');%X轴标签
set(HH1,'LineStyle','-','LineWidth',2,'color','b');%设置左坐标系曲线线型、标记及颜色
set(HH2,'LineStyle','-.','LineWidth',2,'color','r');%设置右坐标系曲线线型、标记及颜色
%以下设置坐标系参数
set(AAX(1),'XColor','k','YColor','b');%设置左坐标系颜色
set(AAX(2),'XColor','k','YColor','r');%设置右坐标系颜色
set(get(AAX(1),'Ylabel'),'string','Curvature Radius\it \rho\rm (mm)','color','b','Fontsize', 14,'fontname','times new roman');%%左侧坐标轴label黑色字号字体
set(get(AAX(2),'Ylabel'),'string','Convex Conditional Function','color','r','Fontsize', 14,'fontname','times new roman');%%右侧坐标轴label字体黑色字号字体
xlabel('Polar Angle of Driving Gear {\itθ} (rad)','color','k','Fontsize', 14,'fontname','times new roman');%设置X轴标签
HH1H=legend([HH1,HH2],{'Curvature Radius','Convex Conditional Function'},'Fontsize', 14);%HHHH1,'Δαc',
%legend('boxoff');  %legend去掉框
%box off;
hold off

figure(4)%从动轮的曲率半径
hold on
plot(sita2,qlbj2,'Color','b','LineStyle','-','linewidth',2);%画从动轮的曲率半径
axis([0 2*pi 0 120]);
ylabel('Curvature Radius\it \rho\rm (mm)','FontName','Times New Roman', 'Fontsize', 14);
set(gca, 'xlim',[0  2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman', 'Fontsize', 14);
set(gca,'YTick',0:20:120)
set(gca,'XTickLabel',{'0','0.5π','π','1.5π','2π'},'Fontsize', 14,'fontname','times new roman');
xlabel('Polar Angle of Driven Gear {\itθ} (rad)','color','k','Fontsize', 14,'fontname','times new roman');%设置X轴标签
% set(gca,'yTickLabel',{'-2000','-1500','-1000','-500','0','500','1000','1500','2000'})
hold off

figure(5)%啮合重合度
hold on
plot(sita1sita1,Chhedu,'Color','b','LineStyle','-','linewidth',2);%画啮合重合度
axis([0 2*pi 1.658   1.663]);
ylabel('Contact Ratio\it \epsilon\rm (mm)','FontName','Times New Roman', 'Fontsize', 14);
xlim([0 2*pi]);
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'YTick',1.658:0.001:1.663)
set(gca,'XTickLabel',{'0','0.5π','π','1.5π','2π'},'Fontsize', 14,'fontname','times new roman');
xlabel('Polar Angle of Driving Gear {\itθ} (rad)','color','k','Fontsize', 14,'fontname','times new roman');%设置X轴标签
% set(gca,'yTickLabel',{'-2000','-1500','-1000','-500','0','500','1000','1500','2000'})



figure(6)%从动轮角速度&从动轮角加速度
set(gca,'Position',[0.17303 0.15455 0.66 0.74]);%设置图幅位置
hold on
[AAAX,HHH1,HHH2] = plotyy(sita2,w2, sita2,alfa2,'plot');
set(AAAX(:),'Ycolor','k','Fontsize', 14); %设定两个Y轴的刻度字号
set(AAAX(1),'ylim',[0.45,0.8],'ytick',0.45:0.05:0.8,'FontName','Times New Roman');  %左Y轴的范围
set(AAAX(2),'ylim',[-0.15,0.2],'ytick',-0.1:0.05:0.2,'FontName','Times New Roman');  %右Y轴的范围
set(AAAX(1),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %左X轴的范围
set(AAAX(2),'xlim',[0,2*pi],'xtick',0:pi/2:2*pi,'FontName','Times New Roman');  %右X轴的范围
set(AAAX,'XTickLabel',{'0','0.5π','π','1.5π','2π'},'Fontsize', 14,'fontname','times new roman');%X轴标签
set(HHH1,'LineStyle','-','LineWidth',2,'color','b');%设置左坐标系曲线线型、标记及颜色
set(HHH2,'LineStyle','-.','LineWidth',2,'color','r');%设置右坐标系曲线线型、标记及颜色
%以下设置坐标系参数
set(AAAX(1),'XColor','k','YColor','b');%设置左坐标系颜色
set(AAAX(2),'XColor','k','YColor','r');%设置右坐标系颜色
set(get(AAAX(1),'Ylabel'),'string','Angular Velocity\it \omega^\prime\rm (rad\cdot s^{-1})','color','b','Fontsize', 14,'fontname','times new roman');%%左侧坐标轴label黑色字号字体
set(get(AAAX(2),'Ylabel'),'string','Angular Acceleration\it \alpha^\prime\rm (rad\cdot s^{-2})','color','r','Fontsize', 14,'fontname','times new roman');%%右侧坐标轴label字体黑色字号字体
xlabel('Polar Angle of Driven Gear {\itθ} (rad)','color','k','Fontsize', 14,'fontname','times new roman');%设置X轴标签
HHH1H=legend([HHH1,HHH2],{'Angular Velocity of Driven Gear','Angular Acceleration of Driven Gear'},'Fontsize', 14);%HHHH1,'Δαc',
%legend('boxoff');  %legend去掉框
%box off;
hold off


figure(7)%主动轮节曲线、齿根圆、齿顶圆
axis equal 
hold on
% plot(r1.*cos(Sitaa1)+20*cos(sita1+u1),r1.*sin(Sitaa1)+20*sin(sita1+u1),'Color','r','LineStyle','-','linewidth',2);%画主动轮节曲线
 plot(r1.*cos(sita1),r1.*sin(sita1),'Color','b','LineStyle','-','linewidth',2);%画主动轮节曲线
 plot(Rf1.*cos(Sitaf1),Rf1.*sin(Sitaf1),'Color','c','LineStyle','-','linewidth',2);%画主动轮齿根圆
 plot(Ra1.*cos(Sitaa1),Ra1.*sin(Sitaa1),'Color','r','LineStyle','-','linewidth',2);%画主动轮齿顶圆
title('主动轮节曲线、齿根圆、齿顶圆') 
hold off
