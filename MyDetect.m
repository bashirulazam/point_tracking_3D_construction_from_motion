function z_upd = MyDetect(s_pred,Sigma_pred,NextFrame,bound_box,L,firstFrame)

Sigma_pred
[V,D] = eig(Sigma_pred(1:2,1:2),'vector');
D = sort(D,'descend');
search_x = round(s_pred(1) - 2*D(1):s_pred(1) + 2*D(1));
search_y = round(s_pred(2) - 2*D(2):s_pred(2) + 2*D(2));
search_x = search_x(search_x > 1 & search_x < 321);
search_y = search_y(search_y > 1 & search_y < 241);
%c = normxcorr2(bound_box,NextFrame(search_y,search_x));
c = normxcorr2(bound_box,NextFrame(search_y,search_x));
c = c(L+1:size(c,1)-L,L+1:size(c,2)-L);
% figure
% imshow(uint8(NextFrame(search_y,search_x)));
[ypeak, xpeak] = find(c==max(c(:)));
 z_upd(2) = search_y(ypeak);
 z_upd(1) = search_x(xpeak);
z_upd = z_upd';
