clc;
clear;
close all;



t = 0.001:0.001:1;
n = numel(t);
N = 10;
M = 11;

x=5*sin(2*pi*80*t);
% x = randn(n);


primary_resp = IMPULSE1(1,[1, 2, 10],0.001,0.1,1)';
secondary_resp = IMPULSE1(1,[1.5, 2.5, 20],0.001,0.1,1)';


subsecondary_resp = secondary_resp;


x_buff = zeros(1,N);
x_af = zeros(1,n);

for i = 1:n
    x_buff = [x(i) x_buff(1:end-1)];
    x_af(i) = sum(primary_resp.*x_buff);
end


x_expended = zeros(1,M);

x_buff2 = zeros(M,N);
y_buff = zeros(M,N);

sub_y_buff = zeros(M,N);
w = zeros(M,N);

mu = 0.259;


for i = 1:n
    for j = 1:M
        if j == 1
            x_expended(j) = x(i);
        end
        if mod(j,2)==0
            x_expended(j) = sin(pi*(j/2)*x(i));
        else
            x_expended(j) = cos(pi*((j-1)/2)*x(i));
        end
    end


    for l = 1:M
        x_buff2(l,:) = [x_expended(l) x_buff2(l,1:end-1)];
        y_buff(l,:) = sum(w(l,:).*x_buff2(l,:));
    end
    y = sum(y_buff);

    y_new(i) = sum(secondary_resp.*y);

    err(i) = x_af(i) - y_new(i);


    for k = 1:M
        sub_y_buff(k,:) = [x_expended(k) sub_y_buff(k,1: end -1)];

        w(k,:) = w(k,:) + mu*(err(i)*(subsecondary_resp.*sub_y_buff(k,:)));
    end


end



plot(x_af, "black");
hold on;
plot(y_new, "yellow");
hold on;
plot(err,'red');
title('FsLMS')
ylabel('Amplitude')
xlabel('Time')
legend("Noise signal", "Counter signal","Error Signal");
hold off;




function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end




