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
#include <chrono>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::endl;
using std::cout;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  cout << "Initializing..." << endl;
  num_particles = 100;  // TODO: Set the number of particles
  
  weights = vector<double>(num_particles,1.0);
  
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
    
  is_initialized = true;
  cout << "Init Complete" << endl;
      

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
  cout << "Predicting..." << endl;
  std::default_random_engine gen;
  std::normal_distribution<double> x_stddev(0.0,std_pos[0]);
  std::normal_distribution<double> y_stddev(0.0,std_pos[1]);
  std::normal_distribution<double> theta_stddev(0.0,std_pos[2]);
  
  if (yaw_rate > 0.00001){
    for (std::vector<Particle>::iterator it = particles.begin() ; it != particles.end(); ++it){

      it->x += velocity/yaw_rate*(sin(it->theta + yaw_rate*delta_t) - sin(it->theta)) + x_stddev(gen);
      it->y += velocity/yaw_rate*(cos(it->theta) - cos(it->theta +yaw_rate*delta_t)) + y_stddev(gen);
      it->theta += yaw_rate*delta_t + theta_stddev(gen);
    }
  }
  else{
    for (std::vector<Particle>::iterator it = particles.begin() ; it != particles.end(); ++it){
      it->x += velocity*(cos(it->theta))*delta_t + x_stddev(gen);
      it->y += velocity*(sin(it->theta))*delta_t + y_stddev(gen);
      it->theta += yaw_rate*delta_t + theta_stddev(gen);
    }
  }
  cout << "Prediction Done" << endl;

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
  if(observations.size() > 0){
    cout << "Updating weights..." << endl;
    for (std::vector<Particle>::iterator it = particles.begin() ; it != particles.end(); ++it){
      it->weight = 1;
      vector<int> ass;
      vector<double> x_ass;
      vector<double> y_ass;
      for(std::vector<LandmarkObs>::const_iterator obs_it = observations.begin(); obs_it != observations.end(); ++obs_it){
        double x_trans = it->x + (cos(it->theta) * obs_it->x) - (sin(it->theta) * obs_it->y);
        x_ass.push_back(x_trans);
        double y_trans = it->y + (sin(it->theta) * obs_it->x) + (cos(it->theta) * obs_it->y);
        y_ass.push_back(y_trans);
        double min_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
        int min_ndx = 0;
        int min_id = map_landmarks.landmark_list[0].id_i;
        for(int i = 1; i < map_landmarks.landmark_list.size(); i++){
          double temp_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f);
          if (temp_dist < min_dist){
            min_dist = temp_dist;
            min_ndx = i;
            min_id = map_landmarks.landmark_list[i].id_i;
          }
        }
        ass.push_back(min_id);
        
        it->weight *= multiv_prob( std_landmark[0], std_landmark[1],x_trans, y_trans, map_landmarks.landmark_list[min_ndx].x_f, map_landmarks.landmark_list[min_ndx].y_f);
        if(it->weight < 0.0){
          cout << it->weight << endl;
        }
        
      }
      SetAssociations(*it, ass, x_ass, y_ass);
      weights[static_cast<int>(it-particles.begin())] = it->weight;

     }


    cout << "Weights updated" << endl;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  cout << "Resampling..." << endl;
  
  vector<Particle> newSample;
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine gen (seed);

  std::uniform_int_distribution<int> rand_ints(0,num_particles);
  
  int ndx = rand_ints(gen);
  cout << "Random_ndx " << ndx << endl;
  double beta = 0.0;
  vector<double>::iterator wMax = max_element(weights.begin(), weights.end());
  for(int i = 0; i < num_particles; i++){
    beta += 2.0 * (*wMax) * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    
    while(weights[ndx] < beta){
      beta -= weights[ndx];
      if(ndx < num_particles - 1){
        ndx++;
      }
      else{
        ndx = 0;
      }
    }
    Particle forPushing = particles[ndx];
    forPushing.id = i;
    newSample.push_back(forPushing);
  }
  particles = newSample;
  cout << "Resample completed" << endl;
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