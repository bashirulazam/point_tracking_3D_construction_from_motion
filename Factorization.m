clear all
close all
clc
%% Factorization Method to reconstruct the 3D structure and motion
% Assumptions: orthographic projection; only one rigid motion

% Known: M image points on each of the N images
N = 30; % 30 images
M = 11; %    image points in one image
load('coords.mat');
pts = coords(1:2, :, :);
%% Calculate the relative 2D coordinates
pts2D = zeros(2, M, N);
for i = 1:N
    [c, r] = relative2D(pts(1,:,i), pts(2,:,i));
    pts2D(1,:,i) = c;
    pts2D(2,:,i) = r;
end
%% Construct the W matrix
W = zeros(2*N, M);
for i = 1:N
    W(2*i-1:2*i, :) = pts2D(:,:,i); % The transpose of the relative coordinates 
end

%% Perform SVD on W
[U, S, V] = svd(W);
U1 = U(:, 1:3);
D1 = S(1:3, 1:3);
V1 = V(:, 1:3);
% The new W
R1 = U1*sqrt(D1);
S1 = (V1*sqrt(D1))';
W1 = R1*S1;

%% Find A = Q*Q' (3*3 matrix)

% Construct R_all A_elements = b
R_all = zeros(3*N, 6);
for i = 1:N
    Ri1 = R1(2*i-1,:);
    Ri2 = R1(2*i,:);
    
    j1 = Ri1(1,1); 
    j2 = Ri1(1,2);
    j3 = Ri1(1,3);
    k1 = Ri2(1,1);
    k2 = Ri2(1,2);
    k3 = Ri2(1,3);
    
    R_all(3*i - 2,:) = [j1^2 2*j1*j2 2*j1*j3 j2^2 2*j2*j3 j3^2];
    R_all(3*i - 1,:) = [k1^2 2*k1*k2 2*k1*k3 k2^2 2*k2*k3 k3^2];
    R_all(3*i, :) = [k1*j1 k1*j2+k2*j1 k1*j3+k3*j1 k2*j2 ...
                    k2*j3+j2*k3 k3*j3];
end 

bi = [1; 1; 0];
b = [];
for i = 1:N
    b = [b;bi];
end

% Solve A
A_elements = R_all\b;

A = zeros(3,3);
A(1,1) = A_elements(1,1);
A(1,2) = A_elements(2,1);
A(1,2) = A(2,1);
A(1,3) = A_elements(3,1);
A(1,3) = A(3,1);
A(2,2) = A_elements(4,1);
A(2,3) = A_elements(5,1);
A(2,3) = A(3,2);
A(3,3) = A_elements(6,1);

%% Perform svd on A to find Q
[Ua, Sa, Va] = svd(A);
Q = Ua * sqrt(Sa);

%% Find R and S uniquely
R = R1 * Q;
S = inv(Q) * S1;

%% Plot the 3D points
x = S(1,:);
y = S(2,:);
z = S(3,:);
scatter3(x,y,z);
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%view(0,-90);
set(gca,'Ydir','reverse')
set(gca,'Zdir','reverse')
line(S(1,1:4), S(2,1:4), S(3,1:4));
line([S(1,1) S(1,6)], [S(2,1) S(2,6)], [S(3,1) S(3,6)]);
line([S(1,2) S(1,5)], [S(2,2) S(2,5)], [S(3,2) S(3,5)]);
line([S(1,4) S(1,7)], [S(2,4) S(2,7)], [S(3,4) S(3,7)]);
line([S(1,3) S(1,8)], [S(2,3) S(2,8)], [S(3,3) S(3,8)]);
line([S(1,8) S(1,9)], [S(2,8) S(2,9)], [S(3,8) S(3,9)]);
line([S(1,9) S(1,10)], [S(2,9) S(2,10)], [S(3,9) S(3,10)]);
line([S(1,11) S(1,10)], [S(2,11) S(2,10)], [S(3,11) S(3,10)]);
line([S(1,11) S(1,5)], [S(2,11) S(2,5)], [S(3,11) S(3,5)]);
line([S(1,6) S(1,5)], [S(2,6) S(2,5)], [S(3,6) S(3,5)]);
line([S(1,4) S(1,10)], [S(2,4) S(2,10)], [S(3,4) S(3,10)]);
line([S(1,1) S(1,8)], [S(2,1) S(2,8)], [S(3,1) S(3,8)]);
line([S(1,7) S(1,10)], [S(2,7) S(2,10)], [S(3,7) S(3,10)]);
line([S(1,7) S(1,8)], [S(2,7) S(2,8)], [S(3,7) S(3,8)]);

%% Recover the rotation matrices
Rotation = zeros(3,3,N);
for i = 1:N
    r1 = R(2*i-1, :);
    r2 = R(2*i, :);
    r3 = cross(r1, r2);
    Rotation(:,:,i) = [r1; r2; r3];
    display(strcat("The rotation matrix of the ",num2str(i),"th image"))
    display([r1; r2; r3])
end

