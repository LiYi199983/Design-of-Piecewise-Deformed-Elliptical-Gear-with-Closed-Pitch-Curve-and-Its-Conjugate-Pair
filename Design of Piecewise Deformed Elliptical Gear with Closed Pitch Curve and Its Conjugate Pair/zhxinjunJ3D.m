
function A1A2=zhxinjunJ3D(k1,XN1,XN2,A1,DN,m11,m12,m13)%�ƽ�ָ�������ľ�%�ƽ�ָ�������ľ�
% syms  sita1
p1=A1*(1-k1^2);
JGE=0.000001;
A11=0;
A12=1000;
i=0;
error=0.00001;
golden=(sqrt(5)-1)/2;
A111 =A12-golden*(A12-A11);       %������ֵ
A112 =A11+golden*(A12-A11);
while A12-A11>error                 %ѭ������
%�����ߵķ������
sita1=0:JGE:2*pi/(XN1*m11*DN);
fun1=p1./(A111-p1-A111*k1*cos(XN1*m11*sita1));
f11=(1/m11)*trapz(fun1)*JGE*m11;%���ֵ�һʽ
sita1= 2*pi/(XN1*m11*DN):JGE:2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN);
fun2=p1./(A111-p1-A111*k1*cos(XN1*m12*(sita1-2*pi/(XN1*m11*DN))+2*pi/DN));
f21=(1/m12)*trapz(fun2)*JGE*m12;%���ֵ�һʽ 
sita1= 2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN):JGE:2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN)+2*pi/(XN1*m13*DN);
fun3=p1./(A111-p1-A111*k1*cos(XN1*m13*(sita1-2*pi/(XN1*m11*DN)-2*pi/(XN1*m12*DN))+2*2*pi/DN));
f31=(1/m13)*trapz(fun3)*JGE*m13;%���ֵ�һʽ
    y1=abs(f11+f21+f31-2*pi/XN2);
%�����ߵķ������
%�����ߵķ������
sita1=0:JGE:2*pi/(XN1*m11*DN);
fun1=p1./(A112-p1-A112*k1*cos(XN1*m11*sita1));
f11=(1/m11)*trapz(fun1)*JGE*m11;%���ֵ�һʽ
sita1= 2*pi/(XN1*m11*DN):JGE:2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN);
fun2=p1./(A112-p1-A112*k1*cos(XN1*m12*(sita1-2*pi/(XN1*m11*DN))+2*pi/DN));
f21=(1/m12)*trapz(fun2)*JGE*m12;%���ֵ�һʽ 
sita1= 2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN):JGE:2*pi/(XN1*m11*DN)+2*pi/(XN1*m12*DN)+2*pi/(XN1*m13*DN);
fun3=p1./(A112-p1-A112*k1*cos(XN1*m13*(sita1-2*pi/(XN1*m11*DN)-2*pi/(XN1*m12*DN))+2*2*pi/DN));
f31=(1/m13)*trapz(fun3)*JGE*m13;%���ֵ�һʽ
%�����ߵķ������
    y2=abs(f11+f21+f31-2*pi/XN2);
        
if y1>y2              %�Ƚϲ����ĺ���ֵ�Ĵ�С   
   A11=A111;                     %���л���
   A111=A112;
   y1=y2;
   A112=A11+golden*(A12-A11);
else  
     A12=A112;
     A112=A111;
     y2=y1;
   A111=A12-golden*(A12-A11);
end
i=i+1;
end       %��������������Ϊֹ��ֹͣ����
A1A2=(A11+A12)/2;
return    %�ƽ�ָ����볤��%�ƽ�ָ����볤��%�ƽ�ָ����볤��


