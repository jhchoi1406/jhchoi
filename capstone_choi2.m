clear ;

T = ones(10,10).*260 ;
T(1,1:end) = 250 ;
T(1:end,1) = 250 ;
T(end,1:end) = 270;
T(2:end,end) = 270;

v = VideoWriter('matice.avi');
open(v);

n = 1;
while n ~= 100
    for i = 2:length(T)-1
        for j = 2:length(T)-1
            T(i,j) = (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))./4 ;
            n = n+1;

            contourf(T)
            colormap(jet)
            axis ij
            axis equal
            colorbar
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
end

close(v);
