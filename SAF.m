clc;
clear;
close all;



t = 0.001:0.001:1;
n = numel(t);
 
N = 10;
E = 11;

x=sin(2*pi*80*t);

% x = randn(n);


primary_resp = IMPULSE1(1,[1, 2, 10],0.001,0.1,1)';
secondary_resp = IMPULSE1(1,[1.5, 2.5, 20],0.001,0.1,1)';


subsecondary_resp = secondary_resp;
path_nonlinearfilter = IMPULSE1(1,[1.5, 2.5, 10],0.001,0.25,1)';


x_buff = zeros(1,N);
x_af = zeros(1,n);

for i = 1:n
    x_buff = [x(i) x_buff(1:end-1)];
    x_af(i) = sum(primary_resp.*x_buff);
end
x_buff2 = zeros(E,N);
sub_y_buff = zeros(E,N);
y_hlp = zeros(1,E);

w = zeros(E,N);
x_prior = zeros(E);


% x_inver = zeros(1,E);
C = 1/6*[-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0];
q = zeros(E,4);

u = 0;
un = [u^3,u^2,u,1];
un_ = zeros(E,4);
un_diff = [3*u^2,2*u,1,0];
un_diff_ = zeros(E,4);


mu = 0.399;

mu_2 = 2.52;

for i = 1:n
    x_prior = [x(i) x_prior(1:end-1)];


    for j = 1:E
        x_buff2(j,:) = [x_prior(j) x_buff2(j,1:end-1)];
        d = sum(w(j,:).*x_buff2(j,:));
        del_x = q(j,2)-q(j,2);
        u = (d/del_x)-floor(d/del_x);
        un_ (j,:) = un;
        un_diff_(j,:) = un_diff;
        y_hlp(j) = (un)*C*q(j,:)';
    end
    y_new(i) = sum(y_hlp);

    err(i) = x_af(i) - y_new(i);

    for k = 1:E
        sub_y_buff(k,:) = [x_prior(k) sub_y_buff(k,1: end -1)];
        w(k,:) = w(k,:) + mu*(err(i)*(subsecondary_resp.*sub_y_buff(k,:)));
        q(k,:) = q(k,:)+(mu_2*err(i)*(C'*(un_(k,:).*path_nonlinearfilter)'))';

    end
    

end






plot(x_af, "black");
hold on;
plot(y_new, "yellow");
hold on;
plot(err, "red");
title('SAF')
ylabel('Amplitude')
xlabel('Time')
legend("Noise signal", "Counter signal","Error signal")
hold off;







function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end




