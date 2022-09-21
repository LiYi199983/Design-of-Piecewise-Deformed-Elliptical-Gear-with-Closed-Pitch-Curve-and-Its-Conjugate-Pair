
clear;clc;close all
z1=33;%�����ֳ���%%%%%%%%�߽���Բ���ֳ���Ϊ�������ĳһ�����ĳ˻���������XN1��������������XN2��Ϊ������
mn=3;%����ģ��
Bc=10*pi/180;%������
pxl1= 0.1;%������ƫ����
mt=mn/cos(Bc);%����ģ��mm
an=20*pi/180;%���߳��ν�
at=atan(tan(an)/cos(Bc));%����ѹ����

XN1=1;XN2=3;%����
DN1=3;      %�ֶ���
m11=1.1       ;%�ɲ�Ϊ����
m12=1.6;            %0.5;�ɲ�Ϊ����
m13=1/(3-1/m11-1/m12);%m11*m12/(3*m11*m12-m11-m12);            %0.5;�ɲ�Ϊ����

tttt=1.5;%%%%%%%%%%%%%%%%%%%%%%%%%%��֤1.5��ӹ����һȦ%%%%%%�ǡ�jiange������
jiange=0.0015; %0.003;%%%%%%%%%%%%%%%%%%%%%%%%%���ʱ��%%%%%��3
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

A1=mt*pi*z1/((m11*temM11*quadl(fun1,0,pi/m11))+(m12*temM12*quadl(fun2,pi/m11,pi/m11+pi/m12))+(m13*temM13*quadl(fun3,pi/m11+pi/m12,2*pi)));% %��ֵ���ַ����������ְ볤��
disp(['A1(�����ְ볤�ᣩ=' num2str(A1)]);

 XN=XN2/XN1;
 A2=A1*sqrt(XN^2-pxl1^2*(XN^2-1));%�Ӷ��ְ볤��
 A1A2=A1+A2;
 disp(['A1A2(���ľࣩ=' num2str(A1A2)]);
 
 z2=L2/(pi*mt);%�Ӷ��ֳ���%�Ӷ��ֳ���%�Ӷ��ֳ���
 disp(['z2(�Ӷ��ֳ�����=' num2str(z2)]);


p1=A1*(1-pxl1^2);
 for i=1:1:tttt/jiange;
 sita1(i)=i*((2*pi)/tttt)*jiange;   
if sita1(i)<=2*pi/(DN1*m11);  
   r1(i)=p1/(1-pxl1*cos(m11*sita1(i)));%����
    r1BX1(i)=r1(i);%����
    FD1=i;
elseif sita1(i)>2*pi/(DN1*m11)&&sita1(i)<=2*pi/(DN1*m12)+2*pi/(DN1*m11);  
  r1(i)=p1/(1-pxl1*cos(m12*(sita1(i)-2*pi/(DN1*m11))+1*2*pi/DN1));%����
    r1BX2(i)=r1(i);%����
    FD2=i;
elseif sita1(i)>2*pi/(DN1*m12)+2*pi/(DN1*m11);  
  r1(i)=p1/(1-pxl1*cos(m13*(sita1(i)-2*pi/(DN1*m11)-2*pi/(DN1*m12))+2*2*pi/DN1));%����
    r1BX3(i)=r1(i);%����
     FD3=i;
end
end
rB1=[r1BX1 r1BX2(fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+1):fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+2*pi/(DN1*m12*jiange)/((2*pi)/tttt)))  r1BX3(fix(2*pi/(DN1*m11*jiange)/((2*pi)/tttt)+2*pi/(DN1*m12*jiange)/((2*pi)/tttt)+1): fix(2*pi/jiange/((2*pi)/tttt)))];





axis equal 
figure(1) %�������ֽ�����
hold on
%polar(sita1,r,'--r');%����������Բ
r1x=rB1.*cos(sita1);%ֱ������X��
r1y=rB1.*sin(sita1);%ֱ������Y��

r1BX1x=rB1(1:FD1).*cos(sita1(1:FD1));
r1BX1y=rB1(1:FD1).*sin(sita1(1:FD1));
plot(r1BX1x,r1BX1y,'Color','r','LineStyle','--','linewidth',2);%%%%%%%%%%%%%%%%%��ֱ��������Բ
r1BX2x(FD1:FD2)=rB1(FD1:FD2).*cos(sita1(FD1:FD2));
r1BX2y(FD1:FD2)=rB1(FD1:FD2).*sin(sita1(FD1:FD2));
plot(r1BX2x(FD1:FD2),r1BX2y(FD1:FD2),'Color','b','linewidth',2);%%%%%%%%%%%%%%%%��ֱ��������Բ
r1BX3x(FD2:FD3)=rB1(FD2:FD3).*cos(sita1(FD2:FD3));
r1BX3y(FD2:FD3)=rB1(FD2:FD3).*sin(sita1(FD2:FD3));
plot(r1BX3x(FD2:FD3),r1BX3y(FD2:FD3),'Color','m','LineStyle','-.','linewidth',2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ֱ��������Բ
r1BXx=[r1BX1x,r1BX2x(FD1:FD2),r1BX3x(FD2:FD3)];%%%%%�����cad
r1BXy=[r1BX1y,r1BX2y(FD1:FD2),r1BX3y(FD2:FD3)];%%%%%�����cad
%mattoacad('mattoacad',rBXx,rBXy); %%%%%�����cad
% line([min(r1x)-20,max(r1x)+20],[0,0],'Color','r','LineStyle','-.','linewidth',1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ֱ��������Բ
% line([0,0],[min(r1y)-20,max(r1y)+20],'Color','r','LineStyle','-.','linewidth',1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ֱ��������Բ
legend({'First segment','Second segment','Third segment'},'Fontsize',16,'fontname','times new roman');
% plot(0,0, 'r+');%�����������ĵ�
ylabel('\ity\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 16);
xlabel('\itx\rm-axle (mm)','FontName','Times New Roman', 'Fontsize', 16);
set(gca,'XTick',-150:15:150,'Fontname', 'Times New Roman', 'Fontsize', 16)
set(gca,'YTick',-300:15:300,'Fontname', 'Times New Roman', 'Fontsize', 16)
hold off

