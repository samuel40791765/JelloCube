/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#define rest_length 0.142857
#define short_shear_rest_length 0.20203030688
#define long_shear_rest_length 0.247435582216867
#define bend_rest_length 0.285714

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */

// hooks law: F = −kHook(|L|−R) L/|L| 
// damping:  F = −kdamp*v
// add hook, damping, and force field together
// border plane equations: x
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
  struct world buffer;
  int i,j,k;
  int res = jello->resolution;
  point force;
  force.x = 0.1;
  force.y = 0.1;
  force.z = -10;

  buffer = *jello; // make a copy of jello

  //go through points
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        point hookforce, dampingforce, collisionforce, forcefieldforce, totalforce;
        hookforce.x = 0; hookforce.y = 0; hookforce.z = 0; 
        dampingforce.x = 0; dampingforce.y = 0; dampingforce.z = 0;
        collisionforce.x = 0; collisionforce.y = 0; collisionforce.z = 0;
        if(res != 0) {
          int forcefield_x = (jello->p[i][j][k].x + 2) / 4 * (res - 1);
          int forcefield_y = (jello->p[i][j][k].y + 2) / 4 * (res - 1);
          int forcefield_z = (jello->p[i][j][k].z + 2) / 4 * (res - 1);
          forcefieldforce = jello->forceField[forcefield_x * res * res + forcefield_y * res + forcefield_z];
          // if(i == 7 && j == 7 && k == 7){
          //   printf("%d %d %d ", forcefield_x, forcefield_y, forcefield_z);
          //   printf("%.2f %.2f %.2f\n", forcefieldforce.x, forcefieldforce.y, forcefieldforce.z);
          // }
        }
        //structural springs
        for(int di=-1; di<2; di++) {
          int ip = i + di;
          point distance;
          if(di == 0 || (ip < 0 || 7 < ip))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[ip][j][k], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - rest_length) * ((distance).z / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[ip][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[ip][j][k].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[ip][j][k].z) * (distance).z) / dist * ((distance).z / dist);
          //abs(jello->p[i][j][k].x - jello->p[ip][j][k].x)
          //printf("%f\n", abs(jello->p[i][j][k].x - jello->p[ip][j][k].x));
        }
        for(int dj=-1; dj<2; dj++) {
          int jp = j + dj;
          point distance;
          if(dj == 0 || (jp < 0 || 7 < jp))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[i][jp][k], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - rest_length) * ((distance).z / dist);
          //dampingforce.x += -1 * jello->dElastic * ((jello->v[i][jp][k].x - jello->v[i][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[i][jp][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[i][jp][k].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[i][jp][k].z) * (distance).z) / dist * ((distance).z / dist);
        }
        for(int dk=-1; dk<2; dk++) {
          int kp = k + dk;
          point distance;
          if(dk == 0 || (kp < 0 || 7 < kp))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[i][j][kp], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - rest_length) * ((distance).z / dist);
          //dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][kp].x - jello->v[i][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[i][j][kp].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[i][j][kp].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[i][j][kp].z) * (distance).z) / dist * ((distance).z / dist);
        }

        // shear springs
        for(int di=-1; di<2; di++) 
          for(int dj=-1; dj<2; dj++)
            for(int dk=-1; dk<2; dk++) 
            {
              int ip = i + di; int jp = j + dj; int kp = k + dk;
              if((di==0 && dj==0) || (dj==0 && dk==0) || (di==0 && dk==0) || (kp < 0 || 7 < kp) || (jp < 0 || 7 < jp) || (ip < 0 || 7 < ip))
                continue;
              point distance;
              pDIFFERENCE(jello->p[i][j][k], jello->p[ip][jp][kp], distance);
              double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
              if(di!=0 && dj!=0 && dk!=0) {
                hookforce.x += -1 * jello->kElastic * (dist - long_shear_rest_length) * ((distance).x / dist);
                hookforce.y += -1 * jello->kElastic * (dist - long_shear_rest_length) * ((distance).y / dist);
                hookforce.z += -1 * jello->kElastic * (dist - long_shear_rest_length) * ((distance).z / dist);
              }
              else {
                hookforce.x += -1 * jello->kElastic * (dist - short_shear_rest_length) * ((distance).x / dist);
                hookforce.y += -1 * jello->kElastic * (dist - short_shear_rest_length) * ((distance).y / dist);
                hookforce.z += -1 * jello->kElastic * (dist - short_shear_rest_length) * ((distance).z / dist);
              }
              //dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][kp].x - jello->v[i][j][k].x) * (distance).x) / dist * ((distance).x / dist);
              dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[ip][jp][kp].x) * (distance).x) / dist * ((distance).x / dist);
              dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[ip][jp][kp].y) * (distance).y) / dist * ((distance).y / dist);
              dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[ip][jp][kp].z) * (distance).z) / dist * ((distance).z / dist);
            } 

        //bend springs
        for(int di=-2; di<3; di+=2) {
          int ip = i + di;
          point distance;
          if(di == 0 || (ip < 0 || 7 < ip))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[ip][j][k], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).z / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[ip][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[ip][j][k].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[ip][j][k].z) * (distance).z) / dist * ((distance).z / dist);
          //abs(jello->p[i][j][k].x - jello->p[ip][j][k].x)
          //printf("%f\n", abs(jello->p[i][j][k].x - jello->p[ip][j][k].x));
        }
        for(int dj=-2; dj<3; dj+=2) {
          int jp = j + dj;
          point distance;
          if(dj == 0 || (jp < 0 || 7 < jp))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[i][jp][k], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).z / dist);
          //dampingforce.x += -1 * jello->dElastic * ((jello->v[i][jp][k].x - jello->v[i][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[i][jp][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[i][jp][k].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[i][jp][k].z) * (distance).z) / dist * ((distance).z / dist);
        }
        for(int dk=-2; dk<3; dk+=2) {
          int kp = k + dk;
          point distance;
          if(dk == 0 || (kp < 0 || 7 < kp))
            continue;
          pDIFFERENCE(jello->p[i][j][k], jello->p[i][j][kp], distance);
          double dist = sqrt((distance).x * (distance).x + (distance).y * (distance).y + (distance).z * (distance).z);
          hookforce.x += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).x / dist);
          hookforce.y += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).y / dist);
          hookforce.z += -1 * jello->kElastic * (dist - bend_rest_length) * ((distance).z / dist);
          //dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][kp].x - jello->v[i][j][k].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.x += -1 * jello->dElastic * ((jello->v[i][j][k].x - jello->v[i][j][kp].x) * (distance).x) / dist * ((distance).x / dist);
          dampingforce.y += -1 * jello->dElastic * ((jello->v[i][j][k].y - jello->v[i][j][kp].y) * (distance).y) / dist * ((distance).y / dist);
          dampingforce.z += -1 * jello->dElastic * ((jello->v[i][j][k].z - jello->v[i][j][kp].z) * (distance).z) / dist * ((distance).z / dist);
        }

        if(1.9999999 < fabs(jello->p[i][j][k].x)) {
          double dist = jello->p[i][j][k].x - 2;
          double abs_dist = fabs(jello->p[i][j][k].x - 2);
          if(jello->p[i][j][k].x < 0) {
            dist = jello->p[i][j][k].x + 2;
            //jello->p[i][j][k].x = -2;
          }
          // else
          //   jello->p[i][j][k].x = 2;
          hookforce.x += -1 * jello->kCollision * (abs_dist) * dist / abs_dist;
          dampingforce.x += -1 * jello->dCollision * ((jello->v[i][j][k].x) * dist) / abs_dist * (dist / abs_dist);
          if((0 < jello->v[i][j][k].x && 0 < jello->p[i][j][k].x) || (0 > jello->v[i][j][k].x && 0 > jello->p[i][j][k].x))
            jello->v[i][j][k].x = 0;
          // if((hookforce.x + dampingforce.x < 0 && jello->p[i][j][k].x < 0) || (hookforce.x + dampingforce.x > 0 && jello->p[i][j][k].x > 0)){
          //   hookforce.x = 0;
          //   dampingforce.x = 0;
          // }
        }
        if(1.9999999 < fabs(jello->p[i][j][k].y)) {
          double dist = jello->p[i][j][k].y - 2;
          double abs_dist = fabs(jello->p[i][j][k].y - 2);
          if(jello->p[i][j][k].y < 0) {
            dist = jello->p[i][j][k].y + 2;
            //jello->p[i][j][k].y = -2;
          }
          // else
          //   jello->p[i][j][k].y = 2;

          hookforce.y += -1 * jello->kCollision * (abs_dist) * dist / abs_dist;
          dampingforce.y += -1 * jello->dCollision * ((jello->v[i][j][k].y) * dist) / abs_dist * (dist / abs_dist);
          if((0 < jello->v[i][j][k].y && 0 < jello->p[i][j][k].y) || (0 > jello->v[i][j][k].y && 0 > jello->p[i][j][k].y))
            jello->v[i][j][k].y = 0;
        }
        if(1.9999999 < fabs(jello->p[i][j][k].z)) {
          double dist = jello->p[i][j][k].z - 2;
          double abs_dist = fabs(jello->p[i][j][k].z - 2);
          if(jello->p[i][j][k].z < 0) {
            dist = jello->p[i][j][k].z + 2;
            //jello->p[i][j][k].z = -2;
          }
          // else
          //   jello->p[i][j][k].z = 2;
          hookforce.z += -1 * jello->kCollision * (abs_dist) * dist / abs_dist;
          dampingforce.z += -1 * jello->dCollision * ((jello->v[i][j][k].z) * dist) / abs_dist * (dist / abs_dist);
          if((0 < jello->v[i][j][k].z && 0 < jello->p[i][j][k].z) || (0 > jello->v[i][j][k].z && 0 > jello->p[i][j][k].z))
            jello->v[i][j][k].z = 0;
        }
        if(jello->planeright){
          if(jello->p[i][j][k].x * jello->a + jello->p[i][j][k].y * jello->b + jello->p[i][j][k].z * jello->c + jello->d <= 0.01) {
            double plane_point_x = ((0 - (0*jello->b + 0*jello->c + jello->d))/jello->a);
            double t = (jello->a*plane_point_x - jello->a*jello->p[i][j][k].x + jello->b*0 - jello->b*jello->p[i][j][k].y + jello->c*0 - jello->c*jello->p[i][j][k].z)
             / (jello->a*jello->a + jello->b*jello->b + jello->c*jello->c);
            double distance = fabs(jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)
             / sqrt(jello->a*jello->a + jello->b*jello->b + jello->c*jello->c);
            // double x_vector = ((0 - (jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d))/jello->a);
            // double y_vector = ((0 - (jello->p[i][j][k].x*jello->a + jello->p[i][j][k].z*jello->c + jello->d))/jello->b);
            // double z_vector = ((0 - (jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->d))/jello->c);
            double x_vector = -1 * t * jello->a;
            double y_vector = -1 * t * jello->b;
            double z_vector = -1 * t * jello->c;
            double total_vector = distance;
            if(x_vector < 0.01) {
              hookforce.x += -1 * jello->kCollision * (fabs(total_vector)) * x_vector / fabs(total_vector);
              dampingforce.x += -1 * jello->dCollision * ((jello->v[i][j][k].x) * x_vector) / fabs(total_vector) * (x_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].x = 0;
                // hookforce.x = 0;
                // dampingforce.x = 0;
                // forcefieldforce.x = 0;
              }
            }
            if(y_vector < 0.01) {
              hookforce.y += -1 * jello->kCollision * (fabs(total_vector)) * y_vector / fabs(total_vector);
              dampingforce.y += -1 * jello->dCollision * ((jello->v[i][j][k].y) * y_vector) / fabs(total_vector) * (y_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].y = 0;
                // hookforce.y = 0;
                // dampingforce.y = 0;
                // forcefieldforce.y = 0;
              }
            }
            if(z_vector < 0.01) {
              hookforce.z += -1 * jello->kCollision * (fabs(total_vector)) * z_vector / fabs(total_vector);
              dampingforce.z += -1 * jello->dCollision * ((jello->v[i][j][k].z) * z_vector) / fabs(total_vector) * (z_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].z = 0;
                // hookforce.z = 0;
                // dampingforce.z = 0;
                // forcefieldforce.y = 0;
              }
            }
          }
        }
        else{
          if(jello->p[i][j][k].x * jello->a + jello->p[i][j][k].y * jello->b + jello->p[i][j][k].z * jello->c + jello->d >= 0.01) {
            double plane_point_x = ((0 - (0*jello->b + 0*jello->c + jello->d))/jello->a);
            double t = (jello->a*plane_point_x - jello->a*jello->p[i][j][k].x + jello->b*0 - jello->b*jello->p[i][j][k].y + jello->c*0 - jello->c*jello->p[i][j][k].z)
             / (jello->a*jello->a + jello->b*jello->b + jello->c*jello->c);
            double distance = fabs(jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)
             / sqrt(jello->a*jello->a + jello->b*jello->b + jello->c*jello->c);
            // double x_vector = ((0 - (jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d))/jello->a);
            // double y_vector = ((0 - (jello->p[i][j][k].x*jello->a + jello->p[i][j][k].z*jello->c + jello->d))/jello->b);
            // double z_vector = ((0 - (jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->d))/jello->c);
            double x_vector = -1 * t * jello->a;
            double y_vector = -1 * t * jello->b;
            double z_vector = -1 * t * jello->c;
            double total_vector = distance;
            if(x_vector > 0.01) {
              hookforce.x += -1 * jello->kCollision * (fabs(total_vector)) * x_vector / fabs(total_vector);
              dampingforce.x += -1 * jello->dCollision * ((jello->v[i][j][k].x) * x_vector) / fabs(total_vector) * (x_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].x = 0;
                // hookforce.x = 0;
                // dampingforce.x = 0;
                // forcefieldforce.x = 0;
              }
            }
            if(y_vector > 0.01) {
              hookforce.y += -1 * jello->kCollision * (fabs(total_vector)) * y_vector / fabs(total_vector);
              dampingforce.y += -1 * jello->dCollision * ((jello->v[i][j][k].y) * y_vector) / fabs(total_vector) * (y_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].y = 0;
                // hookforce.y = 0;
                // dampingforce.y = 0;
                // forcefieldforce.y = 0;
              }
            }
            if(z_vector > 0.01) {
              hookforce.z += -1 * jello->kCollision * (fabs(total_vector)) * z_vector / fabs(total_vector);
              dampingforce.z += -1 * jello->dCollision * ((jello->v[i][j][k].z) * z_vector) / fabs(total_vector) * (z_vector / fabs(total_vector));
              if((0 < jello->v[i][j][k].x && 0.0001 < jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d) || 
                (0 > jello->v[i][j][k].x && -0.0001 > jello->p[i][j][k].x*jello->a + jello->p[i][j][k].y*jello->b + jello->p[i][j][k].z*jello->c + jello->d)){
                jello->v[i][j][k].z = 0;
                // hookforce.z = 0;
                // dampingforce.z = 0;
                // forcefieldforce.y = 0;
              }
            }
          }
        }
        // for(int di=-1; di<3; di++) 
        //   for(int dj=-1; dj<3; dj++)
        //     for(int dk=-1; dk<3; dk++) 
        // {
        //   ip = i + di;
        //   jp = j + dj;
        //   kp = k + dk;
        //   if()
        // }
        // point dampingforce;
        // dampingforce.x = -1 * jello->dElastic * jello->v[i][j][k].x;
        totalforce.x = hookforce.x + dampingforce.x + forcefieldforce.x;
        totalforce.y = hookforce.y + dampingforce.y + forcefieldforce.y;
        totalforce.z = hookforce.z + dampingforce.z + forcefieldforce.z;
        a[i][j][k].x = totalforce.x / jello->mass;
        a[i][j][k].y = totalforce.y / jello->mass;
        a[i][j][k].z = totalforce.z / jello->mass;
        // jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        // jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        // jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        // jello->forceField[i * resolution * resolution + j * resolution + k];
      }

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
        // if(2 < jello->p[i][j][k].x)
        //   jello->p[i][j][k].x = 2;
        // else if(jello->p[i][j][k].x < -2)
        //   jello->p[i][j][k].x = -2;
        // if(2 < jello->p[i][j][k].y)
        //   jello->p[i][j][k].y = 2;
        // else if(jello->p[i][j][k].y < -2)
        //   jello->p[i][j][k].y = -2;
        // if(2 < jello->p[i][j][k].z)
        //   jello->p[i][j][k].z = 2;
        // else if(jello->p[i][j][k].z < -2)
        //   jello->p[i][j][k].z = -2;
      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
        // if(2 < jello->p[i][j][k].x)
        //   jello->p[i][j][k].x = 2;
        // else if(jello->p[i][j][k].x < -2)
        //   jello->p[i][j][k].x = -2;
        // if(2 < jello->p[i][j][k].y)
        //   jello->p[i][j][k].y = 2;
        // else if(jello->p[i][j][k].y < -2)
        //   jello->p[i][j][k].y = -2;
        // if(2 < jello->p[i][j][k].z)
        //   jello->p[i][j][k].z = 2;
        // else if(jello->p[i][j][k].z < -2)
        //   jello->p[i][j][k].z = -2;
      }

  return;  
}
