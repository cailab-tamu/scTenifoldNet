function X = read_coil_images(coil_path, obj_num)
% Read COIL-100 images, given an object number
X = zeros(72,128,128);
for i = 1:72
    degrees = 5*(i-1);
    filename = sprintf('%s/obj%d__%d.png',coil_path,obj_num,degrees);
    X(i,:,:) = rgb2gray(imread(filename, 'png'));
end

end

