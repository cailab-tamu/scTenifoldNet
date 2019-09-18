
% change this path to the correct COIL-100 directory
coil_path = '~/Data/coil-100';


cats = read_coil_images(coil_path,14);
dogs = read_coil_images(coil_path,17);

n = size(cats,1);

dim = size(cats,2)*size(cats,3);
X = zeros(dim,n);
Y = zeros(dim,n);
for i = 1:72
    X(:,i) = reshape(cats(i,:,:),[],1);
    Y(:,i) = reshape(dogs(i,:,:),[],1);
end

train_idx = 1:5:72;
test_idx  = setdiff(1:72,train_idx);
num_train = length(train_idx);
num_test  = n - num_train;

%% manifold warping
[P,rX,rY] = manifold_warping(X,Y,'embed',2);
plot_correspondences(rX,rY,P);
