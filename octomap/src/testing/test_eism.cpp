
#include <stdio.h>
#include <octomap/octomap.h>
#include <octomap/math/Utils.h>
#include <octomap/ColorOcTree.h>
#include "testing.h"

using namespace std;
using namespace octomap;
using namespace octomath;

#define sqrt2pi 2.506628274631
#define P_occ_min 0.0000000001
#define P_occ_max 0.9999999999

double ForwardSensorModel(double x, double mu, double sig, double mu_max){
    double p;

    if(x >= mu_max){
        if(mu >= mu_max){
            p = 1.0;
            printf("This never happens, does it?\n");
        }
        else
            p = 0.0;
    }
    else
        p = 1/(sig*sqrt2pi)*exp(-(x-mu)*(x-mu)/(2*sig*sig));;

    return p;
}

double RayInverseSensorModel(std::vector<double>& map_est, std::vector<double>& ray_depths, std::vector<double>& map_est_ISM_r, double sig, double z){
    int n_r = map_est.size();// & INT_MAX;
//    printf("Inside RISM: n_r = %d\n", n_r);
    double P_rplus[n_r+1];
    double notP = 1.0;
    int k;
    for(k = 0; k < n_r; k++){
        P_rplus[k] = notP*map_est[k];
        notP *= (1-map_est[k]);
    }
    P_rplus[n_r] = notP;




    // Obtain unnormalized probabilities and the normalizer
    double tildeP[n_r], a_temp, inv_eta;
    inv_eta = 0.0;
    for(k = 0; k < n_r; k++){
        a_temp = P_rplus[k]*ForwardSensorModel(z, ray_depths[k], sig, 4.0);
        tildeP[k] = map_est[k]*inv_eta+a_temp;
        inv_eta += a_temp;
    }
    inv_eta += P_rplus[n_r]*ForwardSensorModel(z, ray_depths[n_r], sig, 4.0);

    // Make sure 'inv_eta' is above a minimum threshold and normalize probabilities
//    double min_val_inv_eta = 0.000000000001;
    double P_occ_eval;
    if(inv_eta <= 0.0){
//        printf("Error: inverse eta of ray inverse sensor model\nis %e.\n", inv_eta);
        for(k = 0; k < n_r; k++){
            map_est_ISM_r[k] = map_est[k];
        }
    }
    else{
        for(k = 0; k < n_r; k++){
            P_occ_eval = tildeP[k]/inv_eta;
            // Truncate probabilities close to 0 and 1
            if(     P_occ_eval < P_occ_min)
                map_est_ISM_r[k] = P_occ_min;
            else if(P_occ_eval > P_occ_max)
                map_est_ISM_r[k] = P_occ_max;
            else
                map_est_ISM_r[k] = P_occ_eval;
        }
    }

    return 0;
}

int main(int argc, char** argv) {


  ColorOcTree color_tree(0.1);
  OcTree tree (0.1);

  std::vector<double> map_est;
  std::vector<double> ray_depths;
  std::vector<double> map_est_ISM_r;
  double sig = 0.1;
  double z = 2;

  //  point3d origin (10.01, 10.01, 10.02);
  point3d origin (0.01f, 0.01f, 0.01f);
  point3d point_on_surface (2.0f, 0.01f, 0.01f);

  cout << "Generating sphere at " << origin << " ..." << endl;

  unsigned sphere_beams = 360;
  double angle = 2.0*M_PI/double(sphere_beams);
  Pointcloud p;
  octomap::KeyRay* keyray = new octomap::KeyRay;






  for (unsigned i=0; i<sphere_beams; i++) {
    // for (unsigned j=0; j<sphere_beams/2.0; j++) {
      p.push_back(origin+point_on_surface);
      // point_on_surface.rotate_IP (0,0,angle);
      tree.computeRayKeys(origin, point_on_surface,*keyray);

      // std::cout<<"size of ray: "<<keyray->size()<<std::endl;

      int k = 0;

      for (octomap::KeyRay::iterator it = keyray->begin(); it != keyray->end(); it++) {
					tree.setNodeValue(*it, float(0.5), false); // insert freespace measurement
          // point3d endpoint tree.KeyToCoord(*it);
          // ColorOcTreeNode* n = color_tree.updateNode(point3d(tree.keyToCoord(*it)), true);
            // n->setColor( 255*map_est_ISM_r[k],0, 255*(1-map_est_ISM_r[k]));

          k++;
    }
    point_on_surface.rotate_IP (0,angle,0);
  }



  for (unsigned i=0; i<sphere_beams; i++) {
    // for (unsigned j=0; j<sphere_beams/2.0; j++) {
      p.push_back(origin+point_on_surface);
      // point_on_surface.rotate_IP (0,0,angle);
      tree.computeRayKeys(origin, point_on_surface,*keyray);

      // std::cout<<"size of ray: "<<keyray->size()<<std::endl;
      map_est.resize(keyray->size());
      map_est_ISM_r.resize(keyray->size());
      ray_depths.resize(keyray->size());

      int k = 0;
      for (octomap::KeyRay::iterator it = keyray->begin(); it != keyray->end(); it++) {
          if(tree.search(*it)){
            map_est[k] = tree.search(*it)->getValue();
            // if(i==0){
            //   std::cout<<"before: "<<map_est[k]<<std::endl;
            // }
          }else{
            map_est[k] = 0.01;
          }
           // insert freespace measurement
          ray_depths[k] = 2.0*float(k)/keyray->size();
          k++;
      }


      RayInverseSensorModel(map_est, ray_depths, map_est_ISM_r, sig, z);

      k = 0;
      for (octomap::KeyRay::iterator it = keyray->begin(); it != keyray->end(); it++) {
					tree.setNodeValue(*it, float(map_est_ISM_r[k]), false); // insert freespace measurement
          // if(i==0){
          //   std::cout<<map_est_ISM_r[k]<<std::endl;
          // }

          // point3d endpoint tree.KeyToCoord(*it);
          ColorOcTreeNode* n = color_tree.updateNode(point3d(tree.keyToCoord(*it)), true);
          // if(k == 0){
            // n->setColor(255*map_est_ISM_r[k],0,0);
          // }else{
            n->setColor( 255*map_est_ISM_r[k],0, 255*(1-map_est_ISM_r[k]));

          // }

          k++;
			// }
    }
    point_on_surface.rotate_IP (0,angle,0);
  }
  // tree.insertPointCloud(p, origin);

  cout << "Writing to sphere.bt..." << endl;
  EXPECT_TRUE(tree.writeBinary("sphere.bt"));


  std::cout << "\nWriting to / from file\n===============================\n";
  std::cout << "Writing color tree to " << std::endl;
  // write color tree
  EXPECT_TRUE(color_tree.write("sphere_color.ot"));

  std::cout << "Test successful\n";
  return 0;
}
