function [pdf] = gaussian_smooth(mat,fact)

x1=[1:1:size(mat,1)];
y1=[1:1:size(mat,2)];
[y,x]=meshgrid(y1,x1);

pdf=zeros(size(mat,1),size(mat,2));
C1=1/(2*pi*(fact^2));
C2=2*(fact^2);

for i=1:size(mat,1)
    for j=1:size(mat,2)
        if mat(i,j) > 0
        tmpx=x-i;
        tmpy=y-j;
        [inds] = find((tmpx.^2+tmpy.^2)<=9);
        dist=tmpx(inds).^2 + tmpy(inds).^2;
        f=C1*exp(-dist./C2);
        pdf(i,j)=sum(sum(f.*mat(inds)))/sum(sum(f));
        end
    end
    i;
end

end

