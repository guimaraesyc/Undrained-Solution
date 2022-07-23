%Calculates the temperature profile using Fourier series

function T = Fourier(x)

t_int = 20;                 % °C
t_ext = 24.7;                 % °C
width = 0.03;               % m
time = 50;                 % s
diffusion = 1.18*10^-6;      % m2/s

T = diffusion*time/width^2;

z = x/width;

% error = 1;
% temp = 0;
% count = 0;
% while abs(error)>10^-30
%     count = count + 1;
%     error = 2*sin(count*z)/count;
%     temp = temp + error;
%     %x(count) = count;
%     %y(count) = temp;
%     if count==10^3
%         %plot(x,y)
%         break
%     end
% end

theta = 0;
for i=1:50
    M_j = i*pi;
    B_i = 2*(-1)^i/M_j;
    C_i = exp(-M_j^2*T);
    theta = theta + C_i*B_i*sin(M_j*z);
end
T = t_int + (t_ext - t_int)*(theta + z);
T = -1357.2*x^2 + 6.975*x + 28.504;
end


%  for i=1:50
%  z=i/50;
%  x(i)=0.03*z;
%  T(i)=Fourier(x(i));
%  end
%  plot(x,T,'ro-')


% clear
% clc
% 
% temp_y=1;
% mult=1;
% for k=1:6
%     for j=1:9
%         time = mult*j/100
%         temp_x=1;
%         for i=1:50
%             count = i/50;
%             x(temp_x,temp_y) = 0.03*count;
%             y(temp_x,temp_y) = time;
%             z(temp_x,temp_y) = Fourier(x(temp_x,temp_y),time);
%             temp_x = temp_x + 1;
%         end
%         temp_y = temp_y + 1;
%     end
% mult = mult*10;
% end
% h=surf(x,y,z)
% %set(h,'EdgeColor','none');
% set(gca, 'YScale', 'log')
% %shading interp    % interpolate colors across lines and faces 
% legend;
% colorbar;