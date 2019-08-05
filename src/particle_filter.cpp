/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  
  for(int i = 0; i < num_particles; i++){
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    particles.push_back(p);    
  }
      

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  for (std::vector<Particle>::iterator it = particles.begin() ; it != particles.end(); ++it){
    std::default_random_engine gen;
    std::normal_distribution<double> x_stddev(0.0,std_pos[0]);
  	std::normal_distribution<double> y_stddev(0.0,std_pos[1]);
  	std::normal_distribution<double> theta_stddev(0.0,std_pos[2]);
  
    it->x += velocity/yaw_rate*(sin(it->theta + yaw_rate*delta_t) - sin(it->theta)) + x_stddev(gen);
    it->y += velocity/yaw_rate*(cos(it->theta) - cos(it->theta +yaw_rate*delta_t)) + y_stddev(gen);
  	it->theta += yaw_rate*delta_t + theta_stddev(gen);
  }  

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (std::vector<Particle>::iterator it = particles.begin() ; it != particles.end(); ++it){
    it->weight = 1;
    for(std::vector<LandmarkObs>::const_iterator obs_it = observations.begin(); obs_it != observations.end(); ++it){
      
      double x_trans = it->x + (cos(it->theta) * obs_it->x) - (sin(it->theta) * obs_it->y);

      double y_trans = it->y + (sin(it->theta) * obs_it->x) + (cos(it->theta) * obs_it->y);
   
      double min_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
      int min_ndx = 0;
      for(int i = 1; i < map_landmarks.landmark_list.size(); i++){
        double temp_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f);
        if (temp_dist < min_dist){
          min_dist = temp_dist;
          min_ndx = i;
        }
      }
      it->weight *= multiv_prob(x_trans, y_trans, map_landmarks.landmark_list[min_ndx].x_f, map_landmarks.landmark_list[min_ndx].y_f, std_landmark[0], std_landmark[1]);
    }
    weights[it->id] = it->weight;
  }
  
  
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::default_random_engine gen;
  std::discrete_distribution<int> rand_ndx(0, num_particles-1);
  
  vector<Particle> newSample;
  
  int ndx = rand_ndx(gen);
  double beta = 0;
  vector<double>::iterator wMax = max_element(weights.begin(), weights.end());
  for(int i = 0; i < num_particles; i++){
    beta += 2.0 * *wMax * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    while(weights[ndx] < beta){
      beta -= weights[ndx];
      if(ndx < num_particles - 1){
        ndx++;
      }
      else{
        ndx = 0;
      }
    }
    newSample.push_back(particles[ndx]);
  }
  particles = newSample;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}