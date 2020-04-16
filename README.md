CSCI 520, Assignment 1

Samuel Chiang

================

In my implementation of the jello cube, I first set up the force reactions of the cube springs. There are two parts of the spring force reactions, Hook's law and the damping force. Hook's law of spring force is F=−kHook(|L|−R) * ->L / (|L|). To compute force we need the spring rest length, distorted spring length, and the length force of the three vectors. Since there was the possibility of springs being originally distorted in the world file, I had to measure and hard code all the original lengthes of the structural, bend, and shear springs for the spring force algorithm.  Next, I used the vector distance between the two end points of the spring to compute the length forces on x, y and z. The damping force equation is ->F = -kDamp (->va - ->vb) . ->L/ |L| * ->L / (|L|). To compute the damping force, we need the velocity of both ends of the spring, the distance between them, and the length force on the three vectors. All the information necessary used to compute this are both in the jello cube information and already calculated in the process of computing our hook force. 
After testing and making sure the jello cube reacted naturally to force, I moved on to the collision detection and collision springs. I used a pretty straight forward method in detecting whether the absolute value of each three vectors were larger than two. If any one of the absolute values of the three vectors were larger than two, that meant a collision had occurred at one of the 6 faces of the outer cube. If detected, I used the same force reaction laws used in the jello cube springs (hook and damping). The original rest length of the collision spring was 0 and I used the distance between the face and the collision point as the distorted spring length. The length vector was also quite easy to compute in these cases, as we could just subtract/add 2 to the current collision point's vector that had collided with a outer cube face. One slight error I had to fix was that the cube's collision acceleration was not large enough to stop the cube from colliding into the outer cube face in time. This would result in the cube protruding a bit out of the collision face, which seemed quite unnatural. To fix this, I disallowed the positive accelerations in case of collisions to outer cube faces, and only allowed accelerations in opposite directions.
Implementing the force field was the next part. I used the current x, y, z positions of the jello cube points to compute the corresponding force field force on the point. forcefield = forceField[i * resolution * resolution + j * resolution + k] (with i, j, k as the x,y,z coordinates calculated as a value between 0 and resolution).
Finally, I fixed some small camera bugs and adjusted the lighting and color of my project. One notable thing that was bothering me was that the lighting points were too close to the outer cube faces. This resulted in an unlit/black jello cube face every time it collided with the outer cube faces. I set the lighting a bit further away to allow the jello cube have a more natural lighting while it bounced around the space. 
For my extra credit, I added the collision effects when my cube went into an inclined plane. The implementation was a bit harder than simply colliding into the cube faces, since I couldn't just subtract/add 2 to the collision point's colliding vector. First, I calculated on which side of the inclined plane the cube originally was and detected when a point of the cube was about to cross onto the other side. Then to calculate the collision vector, I calculated both the distance from collision point to the plane and the projection of the point onto the plane. Using the vector of the projection onto the plane, I could successfully calculate the three collision vectors (x,y,z) from the collision point and the plane. Using the collision vectors, the hook and damping collision force could be calculated to achieve a natural bounce back when the jello cube collided with an inclined plane.   
My own created world file spins and shoots the cube towards the x=-2 plane with the cube's preset velocity. I also added a skewed corner on one of the corners of the cube, so when the animation first starts, the cube bounces around on its own due to the force of the corner bouncing back. I also added a rotating force field in my world file, so an invisible force continues to take the cube around the outer cube dimension. The animation jpg files attached show the first 20 seconds of my jello cube world file. 


