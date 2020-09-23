function [s_upd,Sigma_upd] = MyUpdate(s_cur,Sigma_cur,NextFrame,bound_box,L,firstFrame)
[s_pred,Sigma_pred] = MyPredict(s_cur,Sigma_cur);
z_upd = MyDetect(s_pred,Sigma_pred,NextFrame,bound_box,L,firstFrame);

H = [1 0 0 0;
     0 1 0 0];
R = [4 0;
     0 4];
K_upd = Sigma_pred*H'*inv(H*Sigma_pred*H' + R);
s_upd = s_pred + K_upd*(z_upd - H*s_pred);
Sigma_upd = (eye(size(Sigma_cur)) - K_upd*H)*inv(Sigma_pred);