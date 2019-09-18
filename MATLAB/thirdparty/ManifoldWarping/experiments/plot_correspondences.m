function plot_correspondences(X1,X2,P,step)

d = size(X1,1);
if nargin < 3, P = zeros(0,2); end
if nargin < 4, step = 1; end
if size(P,2) ~= 2, error('invalid warping path (P)'); end
if d < 2, error('invalid dimensionality (must be  >= 2)'); end
if d > 3
    X1 = X1(1:3,:);
    X2 = X2(1:3,:);
    d = 3;
end


hold all
for i = 1:step:size(P,1)
    a = X1(:,P(i,1))';
    b = X2(:,P(i,2))';
    if d == 2
        plot([a(1),b(1)],[a(2),b(2)],'g--');
    else
        plot3([a(1),b(1)],[a(2),b(2)],[a(3),b(3)],'g--');
    end
end


if d == 2
    plot(X1(1,:),X1(2,:),'-',X2(1,:),X2(2,:),'-');
else
    plot3(X1(1,:),X1(2,:),X1(3,:),'-',X2(1,:),X2(2,:),X2(3,:),'-');
end
axis square

end
