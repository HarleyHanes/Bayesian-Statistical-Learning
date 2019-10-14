ReadMe- NonsquareMatrixTest Solver

10/14/19 Notes
-Redid Scaling function as ScaleSolve. Am getting good results for error in residuals except for situations    where all other tests give near the exact answer, ScaleSolve has error in residuals on order of 10^(-5).
-I am slightly confused about whether the norm of x should be 1. I haven't been able to get it to be norm 1    but A*x is norm 1 so I think what I have is correct and x isn't neccesarily supposed to be norm 1.
-Next Step: Integrating SVD/ MPP/ Ridge decomp. I think that may help in scenarios where those are accurate but    ScaleSolve is not.

10/7/19 Notes
-Created RidgeScale function that performs Ridge decomposition, scales the 
    equation, solves for x, and then unscales x
-I'm unsure if I'm doing the RidgeScaled correctly. I don't know if I scale
    the matrices first and then perform Ridge decomposition or Ridge decomposition
    then scaling. I've been testing code where I do Ridge decomposition first 
    because I noticed that if you scale A and then add the lambda's the A 
    is no longer scaled correctly. Will ask for clarification in class.
Results: RidgeScaled performs well for singular matrices where Ridge where
         Ridge fails but worse where Ridge performs very well

10/4/19 Notes
-Rewrote ridge decomp with matrices rather than ridge function
-Adjusted MPP so correct for matrices with more rows than columns
Results Notes
Overall notes- MPP sometimes beats SVD but SVD is never horribly off like
MPP, Ridge is never much better than MPP so I think there may be an issue
there
Col > Row
    %Nonsingular- All tests performing well except ridge
    %colSingular- All tests performing well except ridge which is getting
                    %residual norms of about 10^(-5)
    %rowSingular- All tests performing as well as A\b except for xMPP which
    %               has 10x the residual of others (O(.1) vs O(1))
    %row&col Singular- All tests performing as well as A\b except MPP which
    %                  has residuals O(10^17)
 Row < Col
    %Nonsingular- All tests performing as well as A\b
    %colSingular- MPP failing with O(10^15)
    %rowSingular- All tests performing as well as A\b 
    %row&col Singular- MPP sometimes fails, sometimes succeeds



10/2/19 Notes
-Expanded code to assess how different sizings and ranks of A affects SVD 
and MPP Decomps along with comparing the left and right eigenvectors of ATA
and AAT to the U and V matrices of SVD(A). See in code comments for findings.
-Changed Ridge lambda to 10^(-6)
-Solver Results
    Worried There may be an issue with my functions because MPP and Ridge
     are failing for full rank matrices
     MPP works for square matrices
    MPP has huge residuals for row singular matrices if there are more
     columns and huge residuls for column singular matrices if there are more
     rows
     ---It gives a very large x value for the singular row.
    SVD is working where A\b works

