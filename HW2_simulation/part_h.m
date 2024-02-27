sys = tf([0.1 8/15],[1 11/15 61/50 53/180 8/150]);
Ts = 0.7356;
u = @(t) sin(t) + 2*sin(2*t) + 3*sin(3*t) + 4*sin(4*t) + 5*sin(5*t) + 6*sin(6*t);
t = 0:0.7356:50;                          
for i=1:length(t)
    U(i) = u(t(i)); 
end

num = get(sys,'num');
den = get(sys,'den');
delay = get(sys,'iodelay');
Az = cell2mat(den);                  
Az = Az./Az(1);                            
if delay >0
Bz = [zeros(1,delay),cell2mat(num)]./Az(1); 
else
Bz = cell2mat(num)./Az(1);  
end

ind = find(Bz == 0);
test = isempty(ind);
if test == 1
        d = 0;
elseif ind(1) == 1
        d = 1;
elseif ind(1)~= 1 
        d = 0;
end
if length(ind)>1
for i=1:length(ind)-1
    if ind(i+1)-ind(i) == 1
       d = i+1;
    else
        break
    end
end
end
B = Bz(d+1:end);                         
A = Az(2:end);                            
nb=length(B)-1;
na=length(A);
nu = na+nb+1;

Y = lsim(sys,U,t);
Theta = zeros(nu,1);                      % Initial Parameters
P = 10^10 * eye(nu,nu);                    % Initial Covariance Matrix
Phi = zeros(1,nu);                        % Initial phi
for i = 1 : length(U)
[a(:,i),b(:,i),P,Theta,Phi,K(:,i)] = RecursiveLeastSquares(U(1:i),Y(1:i),d,nb,na,P,Theta,Phi,i);
end
%% Plots
% Convergence of parameters
Sa = size(a);
Sb = size(b);
colors = ['b','r','g','k'];
figure
for m = 1:Sa(1)
    subplot(Sa(1), 1, m); % Creating subplots in a column for each 'm' value
    plot(t, a(m, :), 'color', colors(m), 'LineWidth', 1.5);
    hold on;
    plot(t, A(m) * ones(length(t), 1), '--', 'color', 'k');
    grid on;
    title(['Convergence of numerator parameters for a_', num2str(m)]);
    xlabel('time (sec.)');
    ylabel('Parameter');
    legend(['a_', num2str(m)], ['A_', num2str(m)]);
end
figure
hold on

for m = 1:Sb(1)
    subplot(Sb(1), 1, m); 
    plot(t, b(m,:), 'color', colors(m), 'LineWidth', 1.5)
    hold on
    plot(t, B(m) * ones(length(t), 1), '--', 'color', 'k')
    grid on
    title(['Convergence of denominator parameters - b_', num2str(m-1)]) % Update title
    xlabel('time (sec.)')
    ylabel('Parameter')
end

y = lsim(tf(b(:,end)',[1 a(:,end)'],Ts,'iodelay',d,'variable','z^-1'),U,t);
% Compare the actual and estimated systems outputs
figure
grid on 
hold on
plot(t,Y(1:end),'--','LineWidth',1.5)
plot(t,y(1:end),'o','LineWidth',1.5,'color','r')
xlabel('time (sec.)')
ylabel('Amplitude')
title('2^n^d Order Estimated System')
legend('System Output','2^n^d Order Estimated Output') 
function [a,b,P,Theta,phi,K] = RecursiveLeastSquares(U,Y,d,nb,na,P,Theta,phi,n)
nu = na+nb+1;                         
for j = 1:nu
        if j <= na 
            if (n-j)<=0
                phi(n,j) = 0;
            else
                phi(n,j) = -Y(n-j);
            end
        else       
            if (n-d-(j-(na+1)))<=0
                phi(n,j) = 0;
            else
                phi(n,j) = U(n-d-(j-(na+1)));
            end
        end
end   
           
           K = P*phi(n,:)'*inv(1+phi(n,:)*P*phi(n,:)');
           Theta = Theta+K*(Y(n)-phi(n,:)*Theta);
           P = P-P*phi(n,:)'*inv(1+phi(n,:)*P*phi(n,:)')*phi(n,:)*P;
                    % Estimated System Parameters
           a = Theta(1:na);
           b = Theta(na+1:end);
end
