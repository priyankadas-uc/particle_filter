function out=my_particle_filter(in,P)

load('magdata.mat');
X=xyz(:,1);
Y=xyz(:,2);
Z=xyz(:,3);
% surf(X,Y,Z);
% plot(rand(10,1))
randx = in(1);
randy = in(2);
x_P=[];
% old = digits(4);

% randx = -87;
% randy = 37.0000;
x=sort(X);
y=sort(Y);
mag_x=[];
mag_y=[];
n1=size(x);
n2=size(y);
est_x=[];
est_y=[];
for i=2:n1(1)
     if(x(i)>randx)
        mag_x=[x(i);x(i-1)];
        break
     end
end
for i=2:n2(1)
     if(y(i)>randy)
        mag_y=[y(i);y(i-1)];
        break
     end
end
pos1 = [];
pos2 = [];

for i = 1:(size(mag_x,1))
    pos1=[pos1 ; find(X==mag_x(i))];
    pos2=[pos2 ; find(Y==mag_y(i))];
end

[val,pos]=intersect(pos1,pos2);
% mag_x
% mag_y
if isempty(val)
    p=find(X>=randx-0.02&X<=randx+0.02);
    n3=size(p);
    for i=1:n3(1)
        est_x=[est_x;X(p(i))];
    end
    pos_est=find(Y>=randy-0.02&Y<=randy+0.02);
    n4=size(pos_est);
    for i=1:n4(1)
        est_y=[est_y;Y(pos_est(i))];
    end
    if n3(1)<n4(1)
        x_P(:,1)=est_x(1:n3(1));
        x_P(:,2)=est_y(1:n3(1));
        x_P(:,3)=3.14/4;
        N=n3(1);
    else
        x_P(:,1)=est_x(1:n4(1));
        x_P(:,2)=est_y(1:n4(1));
        x_P(:,3)=3.14/4;
        N=n4(1);
    end
else
    q=find(Z==Z(val(1)));
    n2=size(q);
    for i=1:n2(1)
        est_x=[est_x;X(q(i))];
        est_y=[est_y;Y(q(i))];
    end
    x_P(:,1)=est_x;
    x_P(:,2)=est_y;
    x_P(:,3)=3.14/4;
    N=n2(1);
end









x=[in(1);in(2);in(3)];
t=in(4);
% if(t==0)
%     x=[10;10;0];
% end

x_N=[0.1;0.1;0.1];
x_R=0.1;
T=101;
% t=0.3;
% N=100;

% data2=xlsread('particle_filter_sensor_odom.csv','AngularVelocity:100');

%initilize our initial, prior particle distribution as a gaussian around
%the true initial value
V=2;
 %define the variance of the initial esimate
% x_P = []; % define the vector of particles

% make the randomly generated particles from the initial prior gaussian distribution
% for i = 1:N
%     x_P(i,1) = x(1) + sqrt(V) * randn;
%     x_P(i,2) = x(2) + sqrt(V) * randn;
%     x_P(i,3) = x(3) + sqrt(V) * randn;
% end
% x_P(1,:)';
%show the distribution the particles around this initial value of x.
% figure(1)
% clf
% subplot(121)
% plot(1,x_P,'.k','markersize',5)
% xlabel('time step')
% ylabel('flight position')
% subplot(122)
% hist(x_P,100)
% xlabel('flight position')
% ylabel('count')
% pause

% %the functions used by the Quail are:
% % x = 0.5*x + 25*x/(1 + x^2) + 8*cos(1.2*(t-1)) + PROCESS NOISE --> sqrt(x_N)*randn
% % z = x^2/20 + MEASUREMENT NOISE -->  sqrt(x_R)*randn;
% 
% %generate the observations from the randomly selected particles, based upon
% %the given function
% z_out = [x^2 / 20 + sqrt(x_R) * randn];  %the actual output vector for measurement values.
 %the actual output vector for measurement values.
x_est = [x]; % time by time output of the particle filters estimate
x_est_out = []; % the vector of particle filter estimates.
% z_out = [data1(1,100); data2(1,100)]
% 
% 
%
F=[1,0,0; 0,1,0; 0,0,1];

% U=[data1(t; -2.1*exp(-7)]

%     %from the previou time step, update the flight position, and observed
%     %position (i.e. update the Quails position with the non linear function
%     %and update from this position what the chasing ninja's see confounded
%     %by the Quails illusions.
%     U=[data1(t); data2(t)];
u=[0.1; 0]
z=[in(1);in(2);in(3)];
z_out=z;
x=x+[u(1)*cos(u(2));u(1)*sin(u(2));u(2)]*t/111000+sqrt(x_N)*randn;
%     x =;0.5*x + 25*x/(1 + x^2) + 8*cos(1.2*(t-1)) +  sqrt(x_N)*randn;
%     z = U + sqrt(x_R)*randn;

%     %Here, we do the particle filter
for i = 1:N
%         %given the prior set of particle (i.e. randomly generated locations
%         %the quail might be), run each of these particles through the state
%         %update model to make a new set of transitioned particles.
    x_P_update(i,:) = x_P(i,:)'+[u(1)*cos(u(2));u(1)*sin(u(2));u(2)]*t/111000+sqrt(x_N)*randn;
%         %with these new updated particle locations, update the observations
%         %for each of these particles.
    z_update(i,:) = x_P_update(i,:)';
%         %Generate the weights for each of these particles.
%         %The weights are based upon the probability of the given
%         %observation for a particle, GIVEN the actual observation.
%         %That is, if we observe a location z, and we know our observation error is
%         %guassian with variance x_R, then the probability of seeing a given
%         %z centered at that actual measurement is (from the equation of a
%         %gaussian)
    z - z_update(i,:)';
    (z - z_update(i,:)')'*(z - z_update(i,:)');
    -((z - z_update(i,:)')'*(z - z_update(i,:)'))/(2*x_R);
    exp(-((z - z_update(i,:)')'*(z - z_update(i,:)'))/(2*x_R));
    P_w(i) = (1/sqrt(2*pi*x_R)) * exp(-((z - z_update(i,:)')'*(z - z_update(i,:)'))/(2));
        
end
%     
%     % Normalize to form a probability distribution (i.e. sum to 1).
P_w = P_w./sum(P_w);
%     
%     %% Resampling: From this new distribution, now we randomly sample from it to generate our new estimate particles
%     
%     %what this code specifically does is randomly, uniformally, sample from
%     %the cummulative distribution of the probability distribution
%     %generated by the weighted vector P_w.  If you sample randomly over
%     %this distribution, you will select values based upon there statistical
%     %probability, and thus, on average, pick values with the higher weights
%     %(i.e. high probability of being correct given the observation z).
%     %store this new value to the new estimate which will go back into the
%     %next iteration
for i = 1 : N
    x_P(i,:) = x_P_update(find(rand <= cumsum(P_w),1),:);
end
%     
%     %The final estimate is some metric of these final resampling, such as
%     %the mean value or variance
x_est = mean(x_P);
% x_est_out(t,:)=x_est;
out=x_est';
%     
%     %{
%     %the after
%     subplot(133)
%     plot(0,x_P_update,'.k','markersize',5)
%     hold on
%     plot(0,x_P,'.r','markersize',5)
%     plot(0,x_est,'.g','markersize',40)
%     xlabel('fixed time point')
%     title('weight based resampling')
%     pause
%     %}
%     % Save data in arrays for later plotting
%     x_out = [x_out x];
%     z_out = [z_out z];
%     x_est_out = [x_est_out x_est];
% %     


end
    
