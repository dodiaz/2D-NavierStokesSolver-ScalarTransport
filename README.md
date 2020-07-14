# 2D_NavierStokesSolver_ScalarTransport

Objective of the project: create a second-order unsteady Navier-Stokes solver on a semi-complex geometry including scalar transport. Pictured below is an image of the geometry

![](geometry.png)


Some notes:
1. There are four new arrays that I created for \Delta u**, \Delta v**, \Delta u*, and \Delta v* and I called them du_ss, dv_ss, du_s, and dv_s for simplicity. IF you think of a better notation, feel free to change it!
2. I've created a variable called h_old which will store H^{n-1} since we will need that for the nth step.


Two things that I'm not sure about that we'll have to think about:
1. the document mentions "you should adjust the outflow velicty so that its bulk value is equal to U (implying that mass is exactly conserved)" but I'm not sure how to do this. 
2. The document doesn't eplain when the scalar transport step should be implemented in the code but it seems like it should be right after step 3 and should be solved at each time step.
