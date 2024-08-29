clear all;
clc;
%problem = randomsystem(3,10,3);
%load('rose.mat')
problem = robspat();
repeats = 1;
error_rschur = 0;
error_rjd = 0;
error_schur = 0;
options = struct();
options.recursive = 'recursive';
[commonroots,options,A] = macaulaylab(problem, 30, options);
for i=1:repeats

    our_solutions = rjd(A);
    rschur_solutions = rschur(A);
    schur_solutions = fixschur(A);
    [accuracy, ~] = residuals(problem,our_solutions,@evalmon);
    error_rjd = error_rjd + accuracy;
    [accuracy, ~] = residuals(problem,rschur_solutions,@evalmon);
    error_rschur = error_rschur + accuracy;
    [accuracy, ~] = residuals(problem,schur_solutions,@evalmon);
    error_schur = error_schur + accuracy;
end
size(A)
size(A{1})
%size(our_solutions)

fprintf("Our method: %e\n",error_rjd / repeats);
fprintf("RSchur method: %e\n",error_rschur/ repeats);
fprintf("Schur method: %e\n",error_schur/ repeats);
save('robspat.mat','A');

