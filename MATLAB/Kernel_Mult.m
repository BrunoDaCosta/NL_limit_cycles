function Kernel = Kernel_Mult(X1,X2,param)

N=size(X1,1);
M=size(X1,2);
g=1;%param.g;


%eight
% g=1;
% l(1)=0.2;
% l(2)=100;


l(1)=0.2;
l(2)=100;                   
X1_(:,1:4)=X1(:,1:4);%-param.mean(:,1:4));
X2_(:,1:4)=X2(:,1:4);%-param.mean(:,1:4));
X1_(:,3:4)=(X1(:,3:4))./sqrt(X1(:,3).^2+X1(:,4).^2);
X2_(:,3:4)=(X2(:,3:4))./sqrt(X2(:,3).^2+X2(:,4).^2);


dist=(X1_(:,1).'-X2_(:,1)).^2+(X1_(:,2).'-X2_(:,2)).^2;
exp_pos=(-dist./(2.*l(1).^2))';
%kernel_pos = (ones(size(dist))-(dist/(l(1).^2)))'; %traingular
%test1=g*exp(exp_pos);

delta_vel=(X1_(:,3).'-X2_(:,3)).^2+(X1_(:,4).'-X2_(:,4)).^2;
exp_vel=(-delta_vel./(2.*l(2).^2))';
%test2=g*exp(exp_vel);

exp_tot=exp_pos+exp_vel;

% X1_(:,3:4)=(X1(:,3:4))./sqrt(X1(:,3).^2+X1(:,4).^2);
% X2_(:,3:4)=(X2(:,3:4))./sqrt(X2(:,3).^2+X2(:,4).^2);
% delta_vel2=(X1_(:,3).'-X2_(:,3)).^2+(X1_(:,4).'-X2_(:,4)).^2;
%new_exp=new_exp+(-delta_vel2./(2.*l(2).^2))';
% test3=-delta_vel2./(2.*l(2).^2)';
% 
% difference1=test2./test1;
% difference1(isnan(difference1))=0;
% difference2=test3./test1;
% difference2(isnan(difference2))=0;

% 
% X1_(:,1:4)=(X1(:,1:4)-param.mean(:,1:4))./(param.sigma(:,1:4));
% X2_(:,1:4)=(X2(:,1:4)-param.mean(:,1:4))./(param.sigma(:,1:4));
% delta_vel1=(X1_(:,3).'-X2_(:,3)).^2+(X1_(:,4).'-X2_(:,4)).^2;
% test1=(-delta_vel1./(2.*l(2).^2))';
% 
% X1_(:,3:4)=(X1(:,3:4))./sqrt(X1(:,3).^2+X1(:,4).^2);
% X2_(:,3:4)=(X2(:,3:4))./sqrt(X2(:,3).^2+X2(:,4).^2);
% delta_vel2=(X1_(:,3).'-X2_(:,3)).^2+(X1_(:,4).'-X2_(:,4)).^2;
% test2=(-delta_vel2./(2.*l(2).^2))';
% 
% test_ratio1=delta_vel2./delta_vel1;
% test_ratio2=test2./test1;


Kernel=g*exp(exp_tot);


% for i = 1:M
%     if(i<2)
% 
%     elseif(i==2)
% 
%     else
%         
%     end
%     if(i<=2)
%         Kernel_M(i,:,:)=(g*exp(-(X1_(:,i).'-X2_(:,i)).^2./(2.*l(1).^2)))';
% %         test=permute(Kernel_M(i,:,:),[2 3 1]);
%         Kernel=Kernel.*permute(Kernel_M(i,:,:),[2 3 1]);
%     else
%         Kernel_M(i,:,:)=(g*exp(-(X1_(:,i).'-X2_(:,i)).^2./(2.*l(2).^2)))';
% %         test=permute(Kernel_M(i,:,:),[2 3 1]);
%         Kernel=Kernel.*permute(Kernel_M(i,:,:),[2 3 1]);
%     end
% end
