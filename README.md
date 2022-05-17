# Bronchial Tree
Codes for the final project of the subject MATLAB-and-its-applications-on-engineering taught in the Universitat Politènica de Catalunya. This project is a part of a more wide research on tuberculosis lesions on animal models of Macacs by the UPC research group BIOCOM-SC.

## Briew description
Take into account that the project is not finished yet.
  - bronchialTree.m : main code. Right now it computes the brochial branches irrigating the lobes of the right lung.
  - branca.m : this function function represents a branch of the bronchial tree given its initial and final point, the diameter, the orientation, the color      and its transparency.
  - angles.m : this function calculates the orientation for the daugther branches given the volumes irrigated for eahc one, the starting point, the mother, bifurcation and separation vectors.
  - punt_pla.m : this function projects a point into a plane defined by its normal vector and a point.
  - segments.mat : 3D matrix filled with integers from 0 to 20 defining each lobe of the lung (1-->20) and the outside of the lung (0).
  - xarxa.mat : is a matrix obtained from the previous one and some addition data from lung TC. The first three columns of the matrix contain cartesina coordinates for each point of 2, this has been computed with addition data that stablished the ditance in milimeteres between points in 1, in the x, y and z directions. The fourth column of this matrix states numbers 1 or 2 if it correspond to the left or right lung, and the fifth column has numbers from 1 to 20 for each lobe.    


Núria Mercadé Besora
