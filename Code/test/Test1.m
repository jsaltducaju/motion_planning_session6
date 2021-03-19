classdef Test1 < matlab.unittest.TestCase
    %This test evaluates if the solver is able to obtain the same
    %solution as the Matlab's library for computing Dubins paths. This
    %assesment is made by comparing the calculated path lengths. 6
    %different scenarios, composed by different starting and final vehicle
    %poses, are used for the test. 
    
    properties
        
    end
    
    methods (Test)
        function testSol(testCase)
            state_is = [0 0 0;  0 0 pi/2;  0 0 pi/2;  0 -1 pi/4; 0 0 pi/4; 0 0 pi/4];
            state_fs = [2 3 pi;  -2.3 1.7 pi;  2.3 1.7 pi;  2 1.2 0; 2 1.2 0; 1.2 1.7 -pi/2];
            [nr_cases,~] = size(state_is);
            for i = 1:nr_cases
                state_i = state_is(i,:)';
                state_f = state_fs(i,:)';
                actSol = HandIn6_function(state_i, state_f);
                expSol = 0.0;
                testCase.verifyEqual(actSol,expSol, 'AbsTol',0.01);
            end
        end 
        function testInputSize(testCase)
            state_i = 1;
            state_f = [2 3 pi]';
            testCase.verifyError(@()HandIn6_function(state_i, state_f),'HandIn6:InputSize')
            state_i = [0 0 0]';
            state_f = 1;
            testCase.verifyError(@()HandIn6_function(state_i, state_f),'HandIn6:InputSize')      
        end
        function testVariableType(testCase)
            state_i = [0 '0' 0]';
            state_f = [2 3 pi]';
            testCase.verifyError(@()HandIn6_function(state_i, state_f),'HandIn6:InputVariableType')
            state_i = [0 0 0]';
            syms f1 f2 f3
            state_f = [f1 f2 f3]';
            testCase.verifyError(@()HandIn6_function(state_i, state_f),'HandIn6:InputVariableType')      
        end
    end
end

