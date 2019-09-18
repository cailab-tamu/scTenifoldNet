
% class TestAlignment(unittest.TestCase):
%   def setUp(self):
%     self.points1 = numpy.array([[0,0],[2,4],[4,4],[6,9]])
%     self.points2 = numpy.array([[0,0],[0,2],[2,4],[6,6]])
% 
%   def test_TrivialAlignment(self):
%     proj = alignment.TrivialAlignment(self.points1, self.points2)
%     p1,p2 = proj.project(self.points1, self.points2)
%     assert_array_equal(self.points1, p1)
%     assert_array_equal(self.points2, p2)
%     p1,p2 = proj.project(self.points1, self.points2, num_dims=1)
%     assert_array_equal(self.points1[:,:1], p1)
%     assert_array_equal(self.points2[:,:1], p2)
% 
%   # TODO: test the other aligners

a=[0,0; 2,4; 4,4; 6,9];
b=[0,0; 0,2; 2,4; 6,9];


