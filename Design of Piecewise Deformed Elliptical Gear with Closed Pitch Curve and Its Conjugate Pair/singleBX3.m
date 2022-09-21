
clear;clc;close all
z1=33;%主动轮齿数%%%%%%%%高阶椭圆齿轮齿数为其阶数和某一奇数的乘积（必须是XN1的整数倍，否则XN2不为整数）
mn=3;%法向模数
Bc=10*pi/180;%螺旋角
pxl1= 0.1;%主动轮偏心率
mt=mn/cos(Bc);%端面模数mm
an=20*pi/180;%刀具齿形角
at=atan(tan(an)/cos(Bc));%端面压力角

XN1=1;XN2=3;%阶数
DN1=3;      %分段数
m11=1.1       ;%可不为整数
m12=1.6;            %0.5;可不为整数
m13=1/(3-1/m11-1/m12);%m11*m12/(3*m11*m12-m11-m12);            %0.5;可不为整数

tttt=1.5;%%%%%%%%%%%%%%%%%%%%%%%%%%保证1.5秒加工完成一圈%%%%%%是“jiange”倍数
jiange=0.0015; %0.003;%%%%%%%%%%%%%%%%%%%%%%%%%间隔时间%%%%%×3
han=1;
cn=0.25;
hat=han*cos(Bc);
ct=cn*cos(Bc);
ha1=han*mn;
ha2=han*mn;
hf1=(han+cn)*mn;
hf2=(han+cn)*mn;

L1=pi*mt*z1;
L2=L1*XN2/XN1;

temM11=sqrt(1+pxl1^2*(m11^2-1))/m11;
temK11=m11^2*pxl1^2/(1+pxl1^2*(m11^2-1));
temM12=sqrt(1+pxl1^2*(m12^2-1))/m12;
temK12=m12^2*pxl1^2/(1+pxl1^2*(m12^2-1));
temM13=sqrt(1+pxl1^2*(m13^2-1))/m13;
temK13=m13^2*pxl1^2/(1+pxl1^2*(m13^2-1));
syms sitaF
fun1=inline(sqrt(1-temK11*sin(sitaF)^2)); 
fun2=inline(sqrt(1-temK12*sin(sitaF)^2)); 
fun3=inline(sqrt(1-temK13*sin(sitaF)^2)); 

A1=mt*pi*z1/((m11*temM11*quadl(fun1,0,pi/m11))+(m12*temM12*quadl(fun2,pi/m11,pi/m11+pi/m12))+(m13*temM13*quadl(fun3,pi/m11+pi/m12,2*pi)));% %数值积分法计算主动轮半长轴
disp(['A1(主动轮半长轴）=' num2str(A1)]);

 XN=XN2/XN1;
 A2=A1*sqrt(XN^2-pxl1^2*(XN^2-1));%从动轮半长轴
 A1A2=A1+A2;
 disp(['A1A2(中心距）=' num2str(A1A2)]);
 
 z2=L2/(pi*mt);%从动轮齿数%从动轮齿数%从动轮齿数
 disp(['z2(从动轮齿数）=' num2str(z2)]);


p1=A1*(1-pxl1^2);
 for i=1:1:tttt/jiange;
 sita1(i)=i*((2*pi)/tttt)*jiange;   
if sita1(i)<=2*pi/(DN1*m11);  
   r1(i)=p1/(1-pxl1*cos(m11*sita1(i)));%极径
    r1BX1(i)=r1(i);%极径
    FD1=i;
elseif sita1(i)>2*pi/(DN1*m11)&&sita1(i)<=2*pi/(DN1*m12)+2*pi/(DN1*m11);  
  r1(i)=p1/(1-pxl1*cos(m12*(sita1(i)-2*pi/(DN1*m11))+1*2*pi/DN1));%极径
    r1BX2(i)=r1(i);%极径
    FD2=i;
elseif sita1(i)>2*pi/(DN1*m12)+2*pi/(DN1*m11);  
  r1(i)=p1/(1-pxl1*cos(m13*(sita1(i)-2*pi/(DN1*m11)-2*pi/(DN1*m12))+2*2*pi/DN1));%极径
    r1BX3(i)=r1(i);%极径
     FD3=i;
end
end
rB1=[r1BX1 r1BX2(fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+1):fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+2*pi/(DN1*m12*jiange)/((2*pi)/tttt)))  r1BX3(fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+2*pi/(DN1*m12*jiange)/((2*pi)/tttt)+1): fix(2*pi/jiange/((2*pi)/tttt)))];





axis equal 
figure(1) %画主动轮节曲线
hold on
%polar(sita1,r,'--r');%画极坐标椭圆
r1x=rB1.*cos(sita1);%直角坐标X轴
r1y=rB1.*sin(sita1);%直角坐标Y轴

r1BX1x=rB1(1:FD1).*cos(sita1(1:FD1));
r1BX1y=rB1(1:FD1).*sin(sita1(1:FD1));
plot(r1BX1x,r1BX1y,'Color','r','LineStyle','--','linewidth',2);%%%%%%%%%%%%%%%%%画直角坐标椭圆
r1BX2x(FD1:FD2)=rB1(FD1:FD2).*cos(sita1(FD1:FD2));
r1BX2y(FD1:FD2)=rB1(FD1:FD2).*sin(sita1(FD1:FD2));
plot(r1BX2x(FD1:FD2),r1BX2y(FD1:FD2),'Color','b','linewidth',2);%%%%%%%%%%%%%%%%画直角坐标椭圆
r1BX3x(FD2:FD3)=rB1(FD2:FD3).*cos(sita1(FD2:FD3));
r1BX3y(FD2:FD3)=rB1(FD2:FD3).*sin(sita1(FD2:FD3));
plot(r1BX3x(FD2:FD3),r1BX3y(FD2:FD3),'Color','m','LineStyle','-.','linewidth',2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%画直角坐标椭圆
r1BXx=[r1BX1x,r1BX2x(FD1:FD2),r1BX3x(FD2:FD3)];%%%%%输出到cad
r1BXy=[r1BX1y,r1BX2y(FD1:FD2),r1BX3y(FD2:FD3)];%%%%%输出到cad
%mattoacad('mattoacad',rBXx,rBXy); %%%%%输出到cad
% line([min(r1x)-20,max(r1x)+20],[0,0],'Color','r','LineStyle','-.','linewidth',1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%画直角坐标椭圆
% line([0,0],[min(r1y)-20,max(r1y)+20],'Color','r','LineStyle','-.','linewidth',1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%画直角坐标椭圆
legend({'First segment','Second segment','Third segment'},'Fontsize',16,'fontname','times new roman');
% plot(0,0, 'r+');%画主动轮中心点
ylabel('\ity\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 16);
xlabel('\itx\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 16);
set(gca,'XTick',-150:15:150,'Fontname', 'Times New Roman', 'Fontsize', 16)
set(gca,'YTick',-300:15:300,'Fontname', 'Times New Roman', 'Fontsize', 16)
hold off

