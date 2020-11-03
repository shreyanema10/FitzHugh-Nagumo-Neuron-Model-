function [I1, I2] = findIthreshold
% find limit cycle I1 and I2
a=0.5;r=0.1;b=0.1;
v0=0.8;w0=0;
vpts=(-1.5:.05:1.5);
I_th = [];
v1 = 1.26;
v2 = 4.73;
% function f'(v)
fdash = @(v) (2*v*(a+1) - 3*v*v -a);
flag = 0;
% Loop over different values of applied current I
for I = 0.2:0.001:2.0  
    %     f = @(v) (-v.*(v-a).*(v-1)+ I);
    %     g = @(v) (b*v)/r;
    %     
    %     eq = @(v) f(v)-g(v);
    %     
    %     fp = solve(eq);
    p = [-1,a+1,-(a+(b/r)),I];
    fp = roots(p);
    fp = fp(imag(fp)==0);
    if (fdash(fp)-r)>0 && flag ==0
        I1 = I;
        flag = 1;
    end
    if (fdash(fp)-r)<0 && flag ==1
        I2 = I;
        flag = 0;
    end
    
    
end
end
    
