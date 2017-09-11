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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


	num_particles = 50;
	is_initialized = true;




	default_random_engine gen;
	double std_x = std[0]; // meters
    double std_y = std[1]; // meters
    double std_theta = std[2]; //meters

	// This line creates a normal (Gaussian) distribution for x, y and psi.
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);

	Particle particle;

	for (int i = 0; i < num_particles; i++) {
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
	}


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	double std_pos_x = std_pos[0];
	double std_pos_y = std_pos[1];
	double std_theta = std_pos[2];

	double x_new;
	double y_new;
	double theta_new;

	for (int i = 0; i < num_particles; i++) {

		if (fabs(yaw_rate) < .000000001){
			x_new = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			y_new = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			theta_new = particles[i].theta + delta_t * yaw_rate;
		} else {
			x_new = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_new = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			theta_new = particles[i].theta + delta_t * yaw_rate;
		}

		normal_distribution<double> dist_pos_x(x_new, std_pos_x);
		normal_distribution<double> dist_pos_y(y_new, std_pos_y);
		normal_distribution<double> dist_theta(theta_new, std_theta);

		particles[i].x = dist_pos_x(gen);
		particles[i].y = dist_pos_y(gen);
		//particles[i].theta = fmod(dist_theta(gen), 2.0*M_PI); // I hear normalization is bad for this simulator
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	int n_obs = observations.size();
	int n_pred = predicted.size();
	for (int i = 0; i < n_obs; i++){
		double x_obs = observations[i].x;
		double y_obs = observations[i].y;
		double min_dist2 = (predicted[0].x - x_obs) * (predicted[0].x - x_obs) + (predicted[0].y - y_obs) * (predicted[0].y - y_obs);
		observations[i].id = predicted[0].id;

		for (int j = 1; j < n_pred; j++){
			double dist2 = (predicted[j].x - x_obs) * (predicted[j].x - x_obs) + (predicted[j].y - y_obs) * (predicted[j].y - y_obs);
			if (dist2 < min_dist2){
				observations[i].id = predicted[j].id;
				min_dist2 = dist2;
			}
		}

	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
	int n_obs = observations.size();
	std::vector<double> log_weights(num_particles);
	double max_log_weight;
	for (int i = 0; i < num_particles; i++){
		std::vector<LandmarkObs> obs_map_coords(n_obs);
		for (int j = 0; j < n_obs; j++){
			obs_map_coords[j].x = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);

			obs_map_coords[j].y = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);
		}

		std::vector<LandmarkObs> landmarks_in_range;
		double min_dist = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);
		int closest_landmark_index = 0;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++){
			double dist_to_landmark = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
			if (dist_to_landmark < sensor_range){
				LandmarkObs landmark;
				landmark.x = map_landmarks.landmark_list[j].x_f;
				landmark.y = map_landmarks.landmark_list[j].y_f;
				landmark.id = map_landmarks.landmark_list[j].id_i;
				landmarks_in_range.push_back(landmark);
			}
			if (dist_to_landmark < min_dist){
				min_dist = dist_to_landmark;
				int closest_landmark_index = j;
			}
		}
		if (landmarks_in_range.size() == 0){
			LandmarkObs landmark;
			landmark.x = map_landmarks.landmark_list[closest_landmark_index].x_f;
			landmark.y = map_landmarks.landmark_list[closest_landmark_index].y_f;
			landmark.id = map_landmarks.landmark_list[closest_landmark_index].id_i;
			landmarks_in_range.push_back(landmark);
		}

		dataAssociation(landmarks_in_range, obs_map_coords);


		double log_weight = 0;
		for (int j = 0; j < n_obs; j++){
			double diff_x = obs_map_coords[j].x - map_landmarks.landmark_list[obs_map_coords[j].id - 1].x_f;
			double diff_y = obs_map_coords[j].y - map_landmarks.landmark_list[obs_map_coords[j].id - 1].y_f;
			log_weight += -(diff_x*diff_x)/(2*std_landmark[0]*std_landmark[0]);
			log_weight += -(diff_y*diff_y)/(2*std_landmark[1]*std_landmark[1]);
		}
		log_weights[i] = log_weight;
		if (i == 0){
			max_log_weight = log_weight;
		} else if (log_weight > max_log_weight){
			max_log_weight = log_weight;
		}

	}

	for (int i = 0; i < num_particles; i++){
		particles[i].weight = exp(log_weights[i] - max_log_weight);
	}
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particle_list;
	vector<double> weights(num_particles);
	for (int i = 0; i < num_particles; i++){

		weights[i] = particles[i].weight;
	}
	default_random_engine gen;
	discrete_distribution<int> dd(weights.begin(), weights.end());

	for (int i = 0; i < num_particles; i++){

		new_particle_list.push_back(particles[dd(gen)]);
	}
	particles = new_particle_list;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
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
