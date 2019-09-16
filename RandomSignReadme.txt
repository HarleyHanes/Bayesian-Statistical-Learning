RandomSigns Change Notes- Harley 9/15
Combined files from ResidualSignTest (Elliot) and RandomSign (Helen) and 
editted as described below. Used RandomSign as the base to edit from and 
copied over sections from ResidualSignTest.
Note: After edditing source code all functions are giving equal z scores
ResidualSignTest
    -Used the method for getting nPositive, nNegative, and u from ResidualSignTest
     so as not to use a for loop but removed find() since we do not use the 
     indices of the changes and it could be computationally expensive
    -Added a +1 to the definition of u since length(find(diff(sign(r)))) counts
     the number of changes in sign so there is 1 more run than that (if you have
     n change points then you have n+1 different sections). 
Z-score calculation
    -Added a +1 to the equation for mean of u from Random Sign(Identicle 
    to Random Sign) to match equation 1.4.1 on p.14 of Hansen
    -Added a sqrt() to the definition of standard deviation since equation 1.4.1
    on p.14 of Hansen is for the sd^2
    -Added function input for user set z-score
Display
    -Adjusted Display message from ResidualSignTest to include significance
     level the test passes/ fails at.
