clc;
clear;
close all;



% t = 0.001:0.001:1;
A = readtable('run1_bandsaw.xlsx');
t = A{:,1};
x = A{:,2};
format long;

n = numel(t);
 
N = 10;
E = 11;

% x=5*sin(2*pi*80*t);

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
y_buff = zeros(1,N);
sub_y_buff = zeros(E,N);
w = zeros(E,N);
% x_inver = zeros(1,E);
C = 1/6*[-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0];
q(:) = zeros(1,4);

u = 0;
un(:) = [u^3,u^2,u,1];
un_diff(:) = [3*u^2,2*u,1,0];


mu = 0.099;

mu_2 = 20.80;
% subplot(3,1,1);
for i = 1:n
    if i<12
        strt = 1;
    else
        strt = i-10;
    end

    for j = strt:i
        x_buff2(j-strt+1,:) = [x(i-j+1) x_buff2(j-strt+1,1:end-1)];
        y_buff(j-strt+1) = sum(w(j-strt+1,:).*x_buff2(j-strt+1,:));
    end
    d = sum(y_buff);
    


    d2 = sum(d);

    del_x = q(2)-q(1);
    u = (d2/del_x)-floor(d2/del_x);
    

    y_new(i) = (un)*C*q';



    err(i) = x_af(i) - y_new(i);

    for k = strt:i
        sub_y_buff(k+1-strt,:) = [x(i-k+1) sub_y_buff(k+1-strt,1: end -1)];
        w(k+1-strt,:) = w(k+1-strt,:) + mu*(err(i)*(subsecondary_resp.*sub_y_buff(k+1-strt,:)));
    end
    q = q+(mu_2*err(i)*(C'*(un.*path_nonlinearfilter)'))';

  
%     plot(q);
%     hold on;
end
% title('Control point after each Iteration')
% ylabel('Magnitude of element of CP')
% xlabel('Elements of Control point')
% hold off;




% subplot(3,1,2);
plot(x_af, "black");
hold on;
plot(y_new, "yellow");
hold on;
plot(err, "red");
title('SAF (Bandsaw Noise)')
ylabel('Amplitude')
xlabel('Time')
legend("Noise signal", "Counter signal","Error Signal")


% subplot(3,1,3);
% plot(x_af, "black");
% hold on;
% plot(err, "red");
% title('Error signal & Counter Signal')
% ylabel('Amplitude')
% xlabel('Time')
% legend("Noise signal","Error signal")



function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end




