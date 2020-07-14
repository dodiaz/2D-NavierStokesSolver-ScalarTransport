clc
clear variables

step1_mat = load('mat1_x.txt');

dt = 0.01;
Re = 400;
dx = 1 / 20;

%% left to right du_ss on full length of domain
n = 499;
A = zeros(n,n);

for i = 1:n
    A(i,i) = 1 + dt / (Re * dx^2);
end

for i = 1: n-1
    A(i, i + 1) = -dt / (2 * Re * dx^2);
    A(i + 1, i) = -dt / (2 * Re * dx^2);
end

A(n,n) = 1 + dt / (2 * Re * dx^2);

du_ss = zeros(size(step1_mat));



% for middle rows
b = step1_mat(2,2:end);

for i = 2:99
    du_ss(i,2:end) = b/A;
end

% for the bottom row and top row
b = step1_mat(1,2:end);

du_ss(1,2:end) = b/A;
du_ss(end,2:end) = b/A;

%% check with code results

du_ss_code = load('du_ss.txt');
du_s_code = load('du_s.txt');


%% Top to bottom triadiag

du_s = zeros(size(step1_mat));

B = zeros(100,100);

for i = 1:100
    B(i,i) = 1 + dt / (Re * dx^2);
end

for i = 1:99
    B(i, i + 1) = -dt / (2 * Re * dx^2);
    B(i + 1, i) = -dt / (2 * Re * dx^2);
end

B(1,1) = 1 + 3 * dt / (2 * Re * dx^2);
B(100,100) = 1 + 3 * dt / (2 * Re * dx^2);


for i = 2:499
    c = du_ss(:,i)';
    du_s(:,i) = c / B;
end

