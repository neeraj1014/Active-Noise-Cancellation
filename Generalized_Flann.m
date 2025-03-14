clc;
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



x_prior = zeros(1,4);
x_exp = zeros(1,18);
y_buff = zeros(18,N);
w = zeros(18,N);

sub_y_buff = zeros(18,N);
y = zeros(1,n);


mu = 0.502;
for i = 1:n
    x_prior = [x(i) x_prior(1:end-1)];

    for j = 1:4
        if j == 1
            x_exp(j) = x_prior(j);
            x_exp(j+1) = sin(pi*x_prior(j));
            x_exp(j+2) = cos(pi*x_prior(j));
        else
            g = 5*j-6;
            x_exp(g) = x_prior(j);
            x_exp(g+1) = sin(pi*x_prior(j));
            x_exp(g+2) = cos(pi*x_prior(j));
            x_exp(g+3) = x_prior(j)*x_exp(2);
            x_exp(g+4) = x_prior(j)*x_exp(3);
        end
    end
    
    for l = 1:18
        y_buff(l,:) = [x_exp(l) y_buff(l,1:end-1)];
        d = w(l,:).*y_buff(l,:);
        D = sum(secondary_resp.*d);
        y(i) = y(i)+D;
    end

    err(i) = x_af(i) - y(i);



    for k = 1:18
        sub_y_buff(k,:) = [x_exp(k) sub_y_buff(k,1:end-1)];
        w(k,:) = w(k,:)+mu*err(i)*(subsecondary_resp.*sub_y_buff(k,:));
    end

end



plot(x_af,'Black');
hold on;
plot(y,'Yellow');
hold on;
plot(err,'Red');
title("Generalized FLANN");
xlabel("Time");
ylabel("Amplitude");
legend("Noise signal","Counter signal","Error signal")
hold off;






function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end