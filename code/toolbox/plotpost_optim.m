function plotpost_optim(U,d,R)


%
%
% Display the resulting density \(f(x,t)\) for \(t\) from 0 to 1.
%
%

save /tmp/all_pd.mat
Jlist                = s2v(R,'J');
Constr               = s2v(R,'Constr');
MinVal               = s2v(R,'Min');
dim                  = numel(d)-1;
V                    = interp(U);
if(dim==1)
  figure(2);
  surf(V(:,:,2));
  shading interp
  camlight
elseif(dim==2)
  sel                = round(linspace(1,d(3),6));
  figure(2);
  imageplot( mat2cell(V(:,:,sel,3),d(1),d(2),ones(6,1)),'', 2,3);axis equal
end


figure(3);
subplot(3,1,1);
plot(Jlist(min(200,numel(Jlist)):end), '.-'); axis tight;
title('J');
subplot(3,1,2);
plot((Constr(5:end)), '.-'); axis tight;
title('div=0 violation');
subplot(3,1,3);
plot(MinVal(5:end), '.-'); axis tight;
title('>0 violation');
fprintf('Vopt finale %4.8f | Minval %4.8f \n',Jlist(end),MinVal(end));


return
figure(4)
for k = 1:size(V,3),surf(V(:,:,k,3));set(gca,'ZLim',[0,12*1e-3],'View',[-0.5 0]),pause(0.1);end

