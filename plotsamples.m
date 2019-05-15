function plotsamples(T,C,B)
    plot(C(:,1),C(:,2));
    hold on;
    plot(T(:,1),T(:,2),'ro');

    sizt = size(T);
    for i = 1:sizt(1)
        plot(B(i,1:2),B(i,3:4));
    end
    
    xlim([0,max(C(:,1))]);
%     ylim([min(T(:,2))-1.5,max(C(:,2))+0.5])
    xlabel('step size \alpha')
    ylabel('objective function f(\alpha)')
end