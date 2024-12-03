close all; clear all; clc

hw06_worker = hw06();

disp('Test p1, option 1')
tol = 1e-8;
hw_assert(abs ( hw06_worker.p1(@f1_test,0, 1, 1, 1) - 1/4) < tol);
hw_assert(abs ( hw06_worker.p1(@f1_test,0, 1, 4, 1) - 0.328125) < tol);
hw_assert(abs ( hw06_worker.p1(@f2_test,0, 1, 1, 1) - 1/8) < tol);
hw_assert(abs ( hw06_worker.p1(@f2_test,0, 1, 4, 1) - 0.2421875) < tol);
hw_assert(abs ( hw06_worker.p1(@f3_test,0, 1, 1, 1) - 1/16) < tol);
hw_assert(abs ( hw06_worker.p1(@f3_test,0, 1, 4, 1) - 0.189697265625) < tol);
hw_assert(abs ( hw06_worker.p1(@f4_test,0, 1, 1, 1) - 1/32) < tol);
hw_assert(abs ( hw06_worker.p1(@f4_test,0, 1, 4, 1) - 0.15393066406250) < tol);

disp('Test p1, option 2')
hw_assert(abs ( hw06_worker.p1(@f1_test,0, 1, 1, 2) - 1/2) < tol);
hw_assert(abs ( hw06_worker.p1(@f1_test,0, 1, 4, 2) - 0.34375) < tol);
hw_assert(abs ( hw06_worker.p1(@f2_test,0, 1, 1, 2) - 1/2) < tol);
hw_assert(abs ( hw06_worker.p1(@f2_test,0, 1, 4, 2) - 0.265625) < tol);
hw_assert(abs ( hw06_worker.p1(@f3_test,0, 1, 1, 2) - 1/2) < tol);
hw_assert(abs ( hw06_worker.p1(@f3_test,0, 1, 4, 2) -  0.220703125) < tol);
hw_assert(abs ( hw06_worker.p1(@f4_test,0, 1, 1, 2) - 1/2) < tol);
hw_assert(abs ( hw06_worker.p1(@f4_test,0, 1, 4, 2) - 0.1923828125) < tol);

disp('Test p1, option 3')
tol = 1e-4;
hw_assert(abs ( hw06_worker.p1(@f1_test,0, 1, 16, 3) - 1/3) < tol);
hw_assert(abs ( hw06_worker.p1(@f2_test,0, 1, 16, 3) - 1/4) < tol);
hw_assert(abs ( hw06_worker.p1(@f3_test,0, 1, 16, 3) - 1/5) < tol);
hw_assert(abs ( hw06_worker.p1(@f4_test,0, 1, 16, 3) - 1/6) < tol);
hw_assert(abs ( hw06_worker.p1(@f5_test,0, 1, 16, 3) - 1/7) < tol);
hw_assert(abs ( hw06_worker.p1(@f6_test,0, 1, 16, 3) - (1-cos(1))) < tol);
hw_assert(abs ( hw06_worker.p1(@f7_test,0, 1, 16, 3) - (exp(1)-1)) < tol);
hw_assert(abs ( hw06_worker.p1(@f8_test,0, 1, 16, 3) - pi/4) < tol);

disp('Test p3, option 1')
tol=1e-8;
hw_assert(abs ( hw06_worker.p3(@f1_test,0, 1, 4, 1) - 1/3) < tol);
hw_assert(abs ( hw06_worker.p3(@f2_test,0, 1, 4, 1) - 1/4) < tol);
hw_assert(abs ( hw06_worker.p3(@f3_test,0, 1, 4, 1) - 1/5) < tol);
hw_assert(abs ( hw06_worker.p3(@f4_test,0, 1, 4, 1) - 1/6) < tol);
hw_assert(abs ( hw06_worker.p3(@f5_test,0, 1, 4, 1) - 1/7) < tol);
hw_assert(abs ( hw06_worker.p3(@f6_test,0, 1, 4, 1) - (1-cos(1))) < tol);
hw_assert(abs ( hw06_worker.p3(@f7_test,0, 1, 4, 1) - (exp(1)-1)) < tol);
hw_assert(abs ( hw06_worker.p3(@f8_test,0, 1, 4, 1) - pi/4) < tol);
disp('Test p3, option 2')
hw_assert(abs ( hw06_worker.p3(@f1_test,0, 1, 4, 2) - 1/3) < tol);
hw_assert(abs ( hw06_worker.p3(@f2_test,0, 1, 4, 2) - 1/4) < tol);
hw_assert(abs ( hw06_worker.p3(@f3_test,0, 1, 4, 2) - 1/5) < tol);
hw_assert(abs ( hw06_worker.p3(@f4_test,0, 1, 4, 2) - 1/6) < tol);
hw_assert(abs ( hw06_worker.p3(@f5_test,0, 1, 4, 2) - 1/7) < tol);
hw_assert(abs ( hw06_worker.p3(@f6_test,0, 1, 4, 2) - (1-cos(1))) < tol);
hw_assert(abs ( hw06_worker.p3(@f7_test,0, 1, 4, 2) - (exp(1)-1)) < tol);
hw_assert(abs ( hw06_worker.p3(@f8_test,0, 1, 4, 2) - pi/4) < tol);
disp('Test p3, option 3')
hw_assert(abs ( hw06_worker.p3(@f1_test,0, 1, 4, 3) - 1/3) < tol);
hw_assert(abs ( hw06_worker.p3(@f2_test,0, 1, 4, 3) - 1/4) < tol);
hw_assert(abs ( hw06_worker.p3(@f3_test,0, 1, 4, 3) - 1/5) < tol);
hw_assert(abs ( hw06_worker.p3(@f4_test,0, 1, 4, 3) - 1/6) < tol);
hw_assert(abs ( hw06_worker.p3(@f5_test,0, 1, 4, 3) - 1/7) < tol);
hw_assert(abs ( hw06_worker.p3(@f6_test,0, 1, 4, 3) - (1-cos(1))) < tol);
hw_assert(abs ( hw06_worker.p3(@f7_test,0, 1, 4, 3) - (exp(1)-1)) < tol);
hw_assert(abs ( hw06_worker.p3(@f8_test,0, 1, 4, 3) - pi/4) < tol);

disp('Test p4, option 3')
hw_assert( test_p4(hw06_worker.p4()) );

% hw06_worker.test_p5();

function ret = f1_test(x)
    ret = x.^2;
end

function ret = f2_test(x)
    ret = x.^3;
end

function ret = f3_test(x)
    ret = x.^4;
end

function ret = f4_test(x)
    ret = x.^5;
end

function ret = f5_test(x)
    ret = x.^6;
end

function ret = f6_test(x)
    ret = sin(x);
end

function ret = f7_test(x)
    ret = exp(x);
end

function ret = f8_test(x)
    ret = 1./(1 + x.^2);
end

function ret = test_p4(xw)
    x = xw(:, 1);
    w = xw(:, 2);
    ret = true;
    for i = 0:2:10
        ret = ret && abs(dot(x.^i, w) - 2/(i+1)) < 1e-12;
        ret = ret && abs(dot(x.^(i+1), w)) < 1e-12;
    end
end

function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end