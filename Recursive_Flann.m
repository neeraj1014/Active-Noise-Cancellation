clc;
clear;
close all;



t = 0.001:0.001:1;
n = numel(t);
 
N = 10;


x=sin(2*pi*80*t);
primary_resp = IMPULSE1(1,[1, 2, 10],0.001,0.1,1)';
secondary_resp = IMPULSE1(1,[1.5, 2.5, 20],0.001,0.1,1)';


subsecondary_resp = secondary_resp;


x_buff = zeros(1,N);
x_af = zeros(1,n);

for i = 1:n
    x_buff = [x(i) x_buff(1:end-1)];
    x_af(i) = sum(primary_resp.*x_buff);
end





y_buff = zeros(4,N);
sub_y_buff = zeros(4,N);
w = zeros(4,N);

x_exp = zeros(1,4);

y = zeros(1,n);

mu = 1.09;
for i = 1:n
    if i > 1
    t = y(i-1);
    else
        t = 0;
    end
    for j = 1:4
        if j == 1
            y_buff(j,:) = [x(i) y_buff(j,1:end-1)];
            x_exp(j) = x(i);
        elseif j == 2
           y_buff(j,:) = [t  y_buff(j,1:end-1)];
            x_exp(j) = t;
        elseif j == 3
            y_buff(j,:) = [sin(t)  y_buff(j,1:end-1)];
            x_exp(j) = sin(t);
        else
            y_buff(j,:) = [cos(t)  y_buff(j,1:end-1)];
            x_exp(j) = cos(t);
        end

        d = ( y_buff(j,:).*w(j,:));
        D = sum(secondary_resp.*d);
        y(i) = y(i)+D;
    end
    err(i) = x_af(i)-y(i);

    for k = 1:4
        sub_y_buff(k, :) = [x_exp(k) sub_y_buff(k, 1:end-1)];
        w(k,:) = w(k,:)+mu*err(i)*(subsecondary_resp.*sub_y_buff(k,:));
       
    end

end


% subplot(2,1,1);
plot(x_af,'Black');
hold on;
plot(y,'Yellow');
hold on;
plot(err,'Red');
title("Recursive FLANN")
xlabel("Time");
ylabel("Amplitude");
legend("Noise signal","Counter signal","Error signal")
hold off;


% subplot(2,1,2);
% plot(x_af,'Black');
% hold on;
% plot(err,'Red');




function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end

