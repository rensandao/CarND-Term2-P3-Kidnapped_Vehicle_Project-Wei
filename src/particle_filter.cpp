/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;
// As recommended from forum;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    //recommended by some guy in Udatciy forum
    num_particles = 101;
  
    //default_random_engine gen;
    normal_distribution<double> dist_x(x,std[0]);
    normal_distribution<double> dist_y(y,std[1]);
    normal_distribution<double> dist_theta(theta,std[2]);
    
    for(int i = 0; i < num_particles; i++)
    {      
   		Particle  temp_particle;

        temp_particle.x  = dist_x(gen);
		temp_particle.y = dist_y(gen);
		temp_particle.theta = dist_theta(gen);
		
		temp_particle.weight= 1.0;
		temp_particle.id = i;

        particles.push_back(temp_particle);
        weights.push_back(1.0);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    /* For the prediction formula about yaw rate is 0, ref to extended kalman filter 20. Sigma a point prediction assignment 1 in lecture 7 */
    
    /*
       NOTICE HERE
       The default_random_engine shall be out the bracket
       NOTICE HERE
     */
    
    //creat new states
	double pred_x ;
    double pred_y;
    double pred_theta;

    for(int i=0; i<num_particles; ++i){
	    //update x, y and the yaw angle when the yaw rate is not equal to zero
        if(fabs(yaw_rate) >= 0.00001){ 
      
                double x_update = velocity * (sin(particles[i].theta + yaw_rate * delta_t) -sin(particles[i].theta));
                x_update = x_update / yaw_rate;
                double y_update = velocity * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
                y_update = y_update / yaw_rate;
                //calculate new states
                pred_x = particles[i].x + x_update;
                pred_y = particles[i].y + y_update;
                pred_theta = particles[i].theta + yaw_rate * delta_t;
        	
        }else{
            
            for(int i=0; i<num_particles; ++i){
                
                pred_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
                pred_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
                pred_theta = particles[i].theta;
            
            }
        }

        normal_distribution<double> dist_x(pred_x,std_pos[0]);
        normal_distribution<double> dist_y(pred_y,std_pos[1]);
        normal_distribution<double> dist_theta(pred_theta,std_pos[2]);
        
        //add random noise 
        pred_x = dist_x(gen);
        pred_y = dist_y(gen);
        pred_theta = dist_theta(gen);
        
        particles[i].x = pred_x;
        particles[i].y = pred_y;
        particles[i].theta = pred_theta;

    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    for(unsigned int i=0; i<observations.size(); ++i){
        
        //initialize distance ranging from min to max
        double min_pred_obs = numeric_limits<double>::max();

        for(unsigned int j=0; j<predicted.size(); ++j){
            
            //calculate the distance between current landmark and predicted landmarks
            double dist_curr_pred = dist(observations[i].x, observations[i].y,predicted[j].x, predicted[j].y);
            //choose the closest landmark, update
            if(dist_curr_pred < min_pred_obs){
                //update the min dist
                min_pred_obs = dist_curr_pred;

                observations[i].id = predicted[j].id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    for(int i = 0; i < num_particles; i ++){

        //initialize particles' coordinate 
        double paticle_x = particles[i].x;
        double paticle_y = particles[i].y;
        double particle_theta = particles[i].theta;
        
        //create vector for storing predicted landmarks within sensor range
        vector<LandmarkObs> Landmarks_pred;

        for(unsigned int j = 0; j<map_landmarks.landmark_list.size(); ++j)
        {
            float Landmark_x = map_landmarks.landmark_list[j].x_f;
            float Landmark_y = map_landmarks.landmark_list[j].y_f;
            int Landmark_id = map_landmarks.landmark_list[j].id_i;
            
            // within the range of sensor measurement        
            if(dist(paticle_x, paticle_y, Landmark_x, Landmark_y) <= sensor_range)
            {
                Landmarks_pred.push_back(LandmarkObs{Landmark_id, Landmark_x, Landmark_y});
            }
        }
        
        //Transform observations to map coordinate
        vector<LandmarkObs> observations_map;
        for(unsigned int j=0; j<observations.size(); ++j)
        {
                        
            //double x_c = observations[j].x;
            //double y_c = observations[j].y;
            //int id_c = observations[j].id;
            
            double map_x = paticle_x + cos(particle_theta) * observations[j].x - sin(particle_theta) * observations[j].y;
            double map_y = paticle_y + sin(particle_theta) * observations[j].x + cos(particle_theta) * observations[j].y;

            observations_map.push_back(LandmarkObs{observations[j].id, map_x,map_y});
        }
        
        //call dataAssociation for landmarks within range and observation landmarks, collect.
        dataAssociation(Landmarks_pred, observations_map);

        double std_x = std_landmark[0];
        double std_y = std_landmark[1];

        //re-initialize weight
        double particle_w = 1.0;    
        particles[i].weight = particle_w;
        
        for(unsigned int j=0; j<observations_map.size(); ++j)
        {
            double obs_x = observations_map[j].x;
            double obs_y = observations_map[j].y;
            double obs_id = observations_map[j].id;
            
            double obs_pred_x = 0.0, obs_pred_y = 0.0;
            for(unsigned int k = 0; k < Landmarks_pred.size(); k ++)
            {
                if(obs_id == Landmarks_pred[k].id)
                {
                    obs_pred_x = Landmarks_pred[k].x;
                    obs_pred_y = Landmarks_pred[k].y;
                 
                }
            }
                        
            //calculate observation particles' weights with multi-variate Gaussian distribution
            double norm_dis = 1.0/ (2.0 * M_PI * std_x * std_y);
            double weight_obs = norm_dis * exp(-(pow(obs_x-obs_pred_x,2)/(2*pow(std_x,2))+(pow(obs_y-obs_pred_y,2)/(2*pow(std_y,2)))));
            
            // undate every particles's weight
            particle_w *= weight_obs;
            particles[i].weight *= weight_obs;
        }

        // undate every particles's weight
        weights[i] = particle_w;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    //default_random_engine gen;
    discrete_distribution<> dist_index(weights.begin(), weights.end());
    vector<Particle> New_Particles;
    
    for(unsigned int i=0; i<particles.size(); ++i){
        int d = dist_index(gen);
        
        New_Particles.push_back(particles[d]);
    }
    
    particles = New_Particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();
	
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
