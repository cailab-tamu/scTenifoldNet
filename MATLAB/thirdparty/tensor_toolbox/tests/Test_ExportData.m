% Testing importing tensors using import_data
classdef Test_ExportData < matlab.unittest.TestCase
    methods (Test)

        function Full(testCase)
            x = import_data('ktensor_small.ktns');
            export_data(x,'ktensor_small_2.ktns');
            x2 = import_data('ktensor_small_2.ktns');
            testCase.verifyEqual(x, x2);
            if exist('ktensor_small_2.ktns', 'file')==2
                delete('ktensor_small_2.ktns');
            end
        end
                    
    end
end
