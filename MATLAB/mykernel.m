function KMN = mykernel(var,resp,param)
    
    if(size(var,2)~=size(resp,2))
        e1=exp(-(pdist2(var(:,1),resp).^2)./(2*param(1)^2));
        e2=exp(-(pdist2(var(:,2),resp).^2)./(2*param(1)^2));
        e3=exp(-(pdist2(var(:,3),resp).^2)./(2*param(2)^2));
        e4=exp(-(pdist2(var(:,4),resp).^2)./(2*param(2)^2));
        KMN=param(3).*e1.*e2.*e3.*e4;
    else
        e1=exp(-(pdist2(var(:,1),resp(:,1)).^2)./(2*param(1)^2));
        e2=exp(-(pdist2(var(:,2),resp(:,2)).^2)./(2*param(1)^2));
        e3=exp(-(pdist2(var(:,3),resp(:,3)).^2)./(2*param(2)^2));
        e4=exp(-(pdist2(var(:,4),resp(:,4)).^2)./(2*param(2)^2));
        KMN=param(3).*e1.*e2.*e3.*e4;
    end
end