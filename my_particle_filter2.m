

x=[10;0;0];

x_N=[0.1;0.1;0.1];
x_R=0.1;
T=100;
u=[0.1;0.1];
N=100;
z_out=[]
x_out=[]
x_real=[0;0;0]
x_real=x;
time=0:0.01:5;

V=2;  %define the variance of the initial esimate
 
x_P = []; % define the vector of particles

% make the randomly generated particles from the initial prior gaussian distribution
for t=1:120
    
    for i = 1:N
        x_P(i,1) = x(1) + sqrt(V) * randn;
        x_P(i,2) = x(2) + sqrt(V) * randn;
        x_P(i,3) = x(3) + sqrt(V) * randn;
    end
 %the actual output vector for measurement values.
     % time by time output of the particle filters estimate
    x_est = x;
     x_est_out = []; % the vector of particle filter estimates.

    F=[1,0,0; 0,1,0; 0,0,1];
    x_real=x_real+[u(1)*cos(x_real(3));u(1)*sin(x_real(3));u(2)]*time(t)

    u=[1; 0.1];

    

    z=x_real+[u(1)*cos(x_real(3));u(1)*sin(x_real(3));u(2)]*time(t)+sqrt(x_N)*randn;
    z_out(t,:)=[x_real(1),x_real(2),x_real(3)];
    x_out(t,:)=[x(1),x(2),x(3)];
%     x =;0.5*x + 25*x/(1 + x^2) + 8*cos(1.2*(t-1)) +  sqrt(x_N)*randn;
%     z = U + sqrt(x_R)*randn;

%     %Here, we do the particle filter
    for i = 1:N
%         %given the prior set of particle (i.e. randomly generated locations
%         %the quail might be), run each of these particles through the state
%         %update model to make a new set of transitioned particles.
        x_P_update(i,:) = x_P(i,:)'+[u(1)*cos(x_real(3));u(1)*sin(x_real(3));u(2)]*time(t);
%         %with these new updated particle locations, update the observations
%         %for each of these particles.
        z_update(i,:) = x_P_update(i,:)';
%         %Generate the weights for each of these particles.
%         %The weights are based upon the probability of the given
%         %observation for a particle, GIVEN the actual observation.
%         %That is, if we observe a location z, and we know our observation error is
%         %guassian w ith variance x_R, then the probability of seeing a given
%         %z centered at that actual measurement is (from the equation of a
%         %gaussian)
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
    x=[x_est(1);x_est(2);x_est(3)];
   
end
plot(x_out(:,1),x_out(:,2),'r')
hold on
plot(z_out(:,1),z_out(:,2),'b')
figure (1)
