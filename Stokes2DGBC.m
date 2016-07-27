function [ P U V X Y] = Stokes2DGBC(g, M, p0, mu, toGraph )

%pu for M = 4:

%[
% p(1,1) 
% p(2,1)
% p(3,1)
% p(4,1)
% p(1,2)
% p(2,2)
% p(3,2)
% p(4,2)
% p(1,3) 
% p(2,3)
% p(3,3)
% p(4,3)
% p(1,4)  
% p(2,4)
% p(3,4)
% p(4,4)
% u(1,1) 
% u(2,1)
% u(3,1)
% u(4,1)
% u(1,2)
% u(2,2)
% u(3,2)
% u(4,2)
% u(1,3) 
% u(2,3)
% u(3,3)
% u(4,3)
% u(1,4)  
% u(2,4)
% u(3,4)
% u(4,4)
%]

%%%%% v = 0, not using this %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v_xx + v_yy = p_y / mu - g/mu
%v(i-1,j) - 4v_(i,j) + v(i+1,j) + v(i,j+1) + v(i,j-1) = d*(p(i+1,j)-p(i-1,j))/(2*mu) + (d^2)*g / mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% not using this either %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pxx + pyy = 0
%p(i-1,j) - 4p_(i,j) + p(i+1,j) + p(i,j+1) + p(i,j-1) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Equation 1
%-p_x / mu + u_xx + u_yy = 0
%-d*(p(i+1,j)-p(i-1,j))/(2*mu) + u(i-1,j) - 4u_(i,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) = 0

%Equation 2
%u_x = 0
%u(i+1,j)-u(i-1,j) = 0

%Solving A * pu = b

%  g = 0;
%  M = 20;
%  p0 = 0;
%  mu = 2;
 toGraph = 1;

A = zeros(2*M*M,2*M*M);
b = zeros(2*M*M,1);

d = (1-0) / (M-1);

%interior points for equation 1
%There will be M*M Equations on the interior for the first equation, minus the BC's
%This will only consume half of the rows for A and b
for j=2:M-1
    for i=2:M-1
        
        %indices for p
        index = sub2ind([M,M],i,j);
        A(index, sub2ind([M,M], i+1,j)) = -1 * d/(2*mu);
        A(index, sub2ind([M,M], i-1,j)) = 1 * d/(2*mu);
             
        %indices for u
        A(index, M^2 + sub2ind([M,M], i-1,j)) = 1;
        A(index, M^2 + index) = -4;
        A(index, M^2 + sub2ind([M,M], i+1,j)) = 1;
        A(index, M^2 + sub2ind([M,M], i,j+1)) = 1;
        A(index, M^2 + sub2ind([M,M], i,j-1)) = 1;
        
        b(index) = 0;
    end
end

%interior points for equation 2
%There will be M*M more equations here, minus the BC's
%This will consume the second half of the rows for A and b
for j=2:M-1
    for i=2:M-1
        
        index = sub2ind([M,M],i,j);        
        
        %indices for u
        A(M^2 + index, M^2 + sub2ind([M,M], i-1,j)) = -1;
        A(M^2 + index, M^2 + sub2ind([M,M], i+1,j)) = 1; 
        
        b(M^2 + index) = 0;
    end
end

%border conditions, set them for the first 16 equations

%p(1,j) = p0
%p(M, j) = 100
for j = 1:M
    index1 = sub2ind([M,M],1,j);
    indexM = sub2ind([M,M],M,j);
       
    %first 16 equations
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = p0;
    b(indexM) = 100;
    
%     %second 16 equations
%     A(M^2 + index1, index1) = 1;
%     A(M^2 + indexM, indexM) = 1;
%     b(M^2 + index1) = p0;
%     b(M^2 + indexM) = 100;
end

%No slip top and bottom
%u(i,1) = 0
%u(i, M) = 0

for i = 2:M-1
    index1 = sub2ind([M,M], i,1);
    indexM = sub2ind([M,M], i,M);
    
    %first 16 equations
    A(index1, M^2 + index1) = 1;
    A(indexM, M^2 + indexM) = 1;
    
    b(index1) = 0;
    b(indexM) = 0;
    
%     %second 16 equations
%     A(M^2 + index1, M^2 + index1) = 1;
%     A(M^2 + indexM, M^2 + indexM) = 1;
%     
%     b(M^2 + index1) = 0;
%     b(M^2 + indexM) = 0;
end


%border conditions, set them for the second 16 equations
%ux(1,j) = 0
%ux(M, j) = 0
%forward difference formula for left side, backward diff. for right side
for j = 1:M
    index1 = sub2ind([M,M],1,j);
    indexM = sub2ind([M,M],M,j);
    
    A(M^2 + index1, M^2 + index1) = -3;
    A(M^2 + index1, M^2 + sub2ind([M,M],2,j)) = 4;
    A(M^2 + index1, M^2 + sub2ind([M,M],3,j)) = -1;
    
    A(M^2 + indexM, M^2 + indexM) = 3;
    A(M^2 + indexM, M^2 + sub2ind([M,M],M-1,j)) = -4;
    A(M^2 + indexM, M^2 + sub2ind([M,M],M-2,j)) = 1;

    b(M^2 + index1) = 0;
    b(M^2 + indexM) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This creates a checkerboard, don't use it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ux(i,1) = 0
% %ux(i,M) = 0
% %forward difference formula for the bottom, backward diff. for the top
% for i = 2:M-1
%     index1 = sub2ind([M,M],i,1);
%     indexM = sub2ind([M,M],i,M);
%     
%     A(M^2 + index1, M^2 + index1) = -3;
%     A(M^2 + index1, M^2 + sub2ind([M,M],i,2)) = 4;
%     A(M^2 + index1, M^2 + sub2ind([M,M],i,3)) = -1;
%     
%     A(M^2 + indexM, M^2 + indexM) = 3;
%     A(M^2 + indexM, M^2 + sub2ind([M,M],i,M-1)) = -4;
%     A(M^2 + indexM, M^2 + sub2ind([M,M],i,M-2)) = 1;
% 
%     b(M^2 + index1) = 0;
%     b(M^2 + indexM) = 0;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This ALSO creates a checkerboard, don't use it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -p_x / mu + u_xx + u_yy = 0
% -d*(p(i+1,j)-p(i-1,j))/(2*mu) + u(i-1,j) - 4u_(i,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) = 0
% for i = 2:M-1
%     
%     index1 = sub2ind([M,M],i,1);
%     indexM = sub2ind([M,M],i,M);
%     
%     %indices for p
% 
%     A(M^2 + index1, sub2ind([M,M], i+1,j)) = -1 * d/(2*mu);
%     A(M^2 + index1, sub2ind([M,M], i-1,j)) = 1 * d/(2*mu);
%     
%     A(M^2 + indexM, sub2ind([M,M], i+1,j)) = -1 * d/(2*mu);
%     A(M^2 + indexM, sub2ind([M,M], i-1,j)) = 1 * d/(2*mu);
% 
%     %indices for u
%     
%     %uxx for top and bottom
%     A(M^2 + index1, M^2 + sub2ind([M,M], i-1,j)) = 1;
%     A(M^2 + index1, M^2 + index) = -2;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i+1,j)) = 1;
%     
%     A(M^2 + indexM, M^2 + sub2ind([M,M], i-1,j)) = 1;
%     A(M^2 + indexM, M^2 + index) = -2;
%     A(M^2 + indexM, M^2 + sub2ind([M,M], i+1,j)) = 1;
%         
%     %Forward difference for uyy on bottom
%     A(M^2 + index1, M^2 + index) = 2;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,2)) = -5;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,3)) = 4;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,4)) = -1;
%        
%     %Backward difference for uyy on top
%     A(M^2 + index1, M^2 + index) = 2;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,M-1)) = -5;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,M-2)) = 4;
%     A(M^2 + index1, M^2 + sub2ind([M,M], i,M-3)) = -1;
%     
%     b(M^2 + index) = 0;
% end

% %py(i,1) = 0
% %py(i,M) = 0
%forward difference formula for the bottom, backward diff. for the top
for i = 2:M-1
    index1 = sub2ind([M,M],i,1);
    indexM = sub2ind([M,M],i,M);
    
    A(M^2 + index1, index1) = -3;
    A(M^2 + index1, sub2ind([M,M],i,2)) = 4;
    A(M^2 + index1, sub2ind([M,M],i,3)) = -1;
    
    A(M^2 + indexM, indexM) = 3;
    A(M^2 + indexM, sub2ind([M,M],i,M-1)) = -4;
    A(M^2 + indexM, sub2ind([M,M],i,M-2)) = 1;

    b(M^2 + index1) = 0;
    b(M^2 + indexM) = 0;
end

%pu = A\b;
pu = pinv(A) * b;
%rcond(A) OR use cond

%transform p and u into 2x2 matrices and graph
P = zeros(M,M);
U = zeros(M,M);
V = zeros(M,M);
X = linspace(0,1,M);
Y = linspace(0,1,M);

for i = 1: M^2
    [tempi tempj] = ind2sub([M,M], i);
    P(tempi, tempj) = pu(i);
    U(tempi, tempj) = pu(i + M^2);
end
    
[X Y] = meshgrid(X,Y);

if(toGraph)
    figure(1);
    surf(Y,X,P);
    %title('Pressure');
    zlabel('p');
    xlabel('x');
    ylabel('y');
    
    figure(2);
    surf(Y,X,U);
    %title('Horizontal Velocity');
    zlabel('u');
    xlabel('x');
    ylabel('y');
    
    figure(3);
    surf(Y,X,V);
    %title('Vertical Velocity');
    zlabel('v');
    xlabel('x');
    ylabel('y');
end

end

