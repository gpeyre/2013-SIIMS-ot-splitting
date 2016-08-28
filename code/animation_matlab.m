close all;
fin=32;



zlevels=0.0005:0.0005:0.016;

for i=1:fin+1,
    
    data=U.M{3}(:,:,i);

    for j=1:32,
        data2(j,:)=data(33-j,:);
    end

    [c2,hc2] = contour(data2,zlevels);
    colormap('jet');
    caxis([0, 0.0173])
    axis square
    set(gca,'XTickLabel','','YTickLabel','','XTick', [],'YTick', [])
    %axis([1 100 1 100]);
 

    if i<10
        
        saveas(gcf, ['anim_' num2str(alpha*100) '_iso_0' num2str(i) '.eps'],'epsc');
    else     
       
        saveas(gcf, ['anim_' num2str(alpha*100) '_iso_' num2str(i) '.eps'],'epsc');
    end
    
    
end

