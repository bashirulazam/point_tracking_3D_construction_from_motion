function [s_pred,Sigma_pred] = MyPredict(s_cur,Sigma_cur);
phi = [1 0 1 0;
       0 1 0 1;
       0 0 1 0;
       0 0 0 1];

Q = [16 0 0 0;
     0 16 0 0;
     0  0 4 0;
     0  0 0 4];
s_pred = phi*s_cur + mvnrnd(zeros(size(s_cur)),sqrt(Q),1)';
Sigma_pred = phi*Sigma_cur*phi' + Q;
