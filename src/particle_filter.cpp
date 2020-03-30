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
using std::normal_distribution;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x,(std[0]));
  normal_distribution<double> dist_y(y,(std[1]));
  normal_distribution<double> dist_theta(theta,(std[2]*1));
  
  for(unsigned int i=0;i<num_particles;i++){
    Particle a;
    a.id = i;
    a.x = dist_x(gen);
    a.y = dist_y(gen);
    a.theta = dist_theta(gen);    
    a.weight = 1.0;
    ParticleFilter::particles.push_back(a);
    
  }
  is_initialized = true;
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
    
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0,std_pos[2]);

  unsigned int i;
  if(yaw_rate == 0){
    for(i=0;i<num_particles;i++){
      ParticleFilter::particles[i].x += ((velocity)*(delta_t)*cos(ParticleFilter::particles[i].theta)) + dist_x(gen);
      ParticleFilter::particles[i].y += (velocity*delta_t*sin(ParticleFilter::particles[i].theta)) + dist_y(gen);
      ParticleFilter::particles[i].theta += dist_theta(gen);
    }
  }else{
    for(i=0;i<num_particles;i++){
      ParticleFilter::particles[i].x += (velocity/yaw_rate)*(sin(ParticleFilter::particles[i].theta + yaw_rate*delta_t)-sin(ParticleFilter::particles[i].theta)) + dist_x(gen);
      ParticleFilter::particles[i].y += (velocity/yaw_rate)*(-cos(ParticleFilter::particles[i].theta + yaw_rate*delta_t)+cos(ParticleFilter::particles[i].theta)) + dist_y(gen);
      ParticleFilter::particles[i].theta += yaw_rate*delta_t + dist_theta(gen);
    }
  }
}

void ParticleFilter::dataAssociation(Particle& particle, vector<LandmarkObs> predicted, 
                                     const vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  
  double new_nearest; // to update the nearest landmark's dist
  int nearestindex; //nearest landmark's id (i.e index i +1)
  vector<int> associat = {}; //associations vector
  vector<double> obs1x={}; //sense_x vector
  vector<double> obs1y={}; //sense_y vector
  unsigned int i,j;
  for(j = 0;j<observations.size();j++){    
    obs1x.push_back(particle.x + (cos(particle.theta) * observations[j].x) - (sin(particle.theta) * observations[j].y));
    obs1y.push_back(particle.y + (sin(particle.theta) * observations[j].x) + (cos(particle.theta) * observations[j].y));
    double nearest_dist = 10000; //setting a arbitrary high value to start the loop
    for(i = 0;i<predicted.size();i++){
      diff = dist(obs1x[j],obs1y[j],predicted[i].x,predicted[i].y);
      if(diff<difflast){
        difflast = diff;
        nearestindex = predicted[i].id;
      }
    }    
    associat.push_back(nearestindex);
   }
  
   ParticleFilter::SetAssociations(particle,associat,obs1x,obs1y);
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

  weights = {}; //resetting weights vector
  double weightsum = 0;  
  for(unsigned int j = 0; j<num_particles;j++){
    vector<LandmarkObs> predicted = {}; //reseting predicted vector

    //creating vector of landmarks within sensor range
    for(unsigned int i=0;i<map_landmarks.landmark_list.size();i++){
      if(dist(ParticleFilter::particles[j].x,ParticleFilter::particles[j].y,map_landmarks.landmark_list[i].x_f,map_landmarks.landmark_list[i].y_f)<sensor_range){
        LandmarkObs a;
        a.id = map_landmarks.landmark_list[i].id_i;
        a.x = map_landmarks.landmark_list[i].x_f;
        a.y = map_landmarks.landmark_list[i].y_f;
        predicted.push_back(a);
      }
    }

    ParticleFilter::dataAssociation(ParticleFilter::particles[j],predicted,observations);
    double b = 1.0;
    ParticleFilter::particles[j].weight = 1.0;
    for(unsigned int k = 0; k<observations.size();k++){
 
      b=multiv_prob(std_landmark[0],std_landmark[1],ParticleFilter::particles[j].sense_x[k],ParticleFilter::particles[j].sense_y[k],map_landmarks.landmark_list[ParticleFilter::particles[j].associations[k]-1].x_f,map_landmarks.landmark_list[ParticleFilter::particles[j].associations[k]-1].y_f);
      
      ParticleFilter::particles[j].weight *= b;
    }
    weights.push_back(ParticleFilter::particles[j].weight);
    weightsum += weights[j];
  }
  //normalizing weights to 0-1
  for(unsigned int m=0;m<weights.size();m++){
    weights[m] = weights[m]/weightsum;
    ParticleFilter::particles[m].weight = weights[m];

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
  uniform_int_distribution<int> dist_x(0,(num_particles-1));
  uniform_int_distribution<int> beta_add(0,100);
  int index = dist_x(gen);
  double beta = 0.0;
  double mw = *max_element(weights.begin(),weights.end());
  for(unsigned int j= 0; j<num_particles;j++){
    beta += beta_add(gen) * 0.01 * 2.0 * mw;
    while(beta > weights[index]){
        beta -= weights[index];
        index = (index + 1) % (num_particles);
      }
    ParticleFilter::particles[j] = ParticleFilter::particles[index];
  }
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