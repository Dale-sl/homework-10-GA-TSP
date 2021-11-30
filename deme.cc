/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include "algorithm"
#include <vector>
#include <iostream>
#include "cities.hh"

using namespace std;

bool compare_fitness(Chromosome*a, Chromosome*b){     //helper function outside of class
  if (a->get_fitness() > b->get_fitness()){
    return true;
  }
  else{
    return false;
  }
} 
// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate)
{
  for(unsigned int i=0; i<pop_size; i++){
    Chromosome* a_chromosome = &Chromosome(cities_ptr);   //make pop_size chromosomes 
    pop_.push_back(a_chromosome);                     //push it back to the population vector 
  }
  mut_rate_ = mut_rate;                             //make the mut_rate_ in the header file equal to the thing in the constructor
}

// Clean up as necessary
Deme::~Deme()
{
  for (int i=0; i<pop_.size(); i++){
    delete(pop_[i]);                    //dealocated everything in the vector and then delete the vector
  }
  //delete(pop_);
}


// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation()
{
  std::vector<Chromosome*> new_children;    //vector constructed for the new children later
  for (int i=0; i<(pop_.size()/2); i++){    //for the population vector with size pop_.size()/2
    Chromosome* parent_x = select_parent();   //get a parent x with the select_parent method
    Chromosome* parent_y = select_parent();   //get a perent y with the select_parent method
    double random_number_x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);   //get a number between 0 and 1 for x

    double random_number_y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);   //get a numher between 0 and 1 for y
    if (mut_rate_ > random_number_x){   //if the mutation rate is greater than the random number selected for x, mutate parent x
      parent_x->mutate();
    }
    if (mut_rate_ > random_number_y){   // if the mutation rate is graeter than the random number selected for y, mutate parent y
      parent_y->mutate();
    }
    std::pair<Chromosome*, Chromosome*> a = parent_x->recombine(parent_y);  //recombine parent x and y to make a pair of two children stored in a
    new_children.push_back(a.first);    //put the first child of the pair in the new vector 
    new_children.push_back(a.second);   //put the second child of the pair in the new vector
    for (int i = 0; i<pop_.size(); i++){  //delete everything in the original population vector and then delete the vector entirely
      delete pop_[i];
    }
    //delete pop_;
    pop_ = new_children;      //make the original population vector of the new children
  }
}

// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const    
{
  Chromosome* c;          
  double max = 0.0; 
  Chromosome* cmax;
  for (int i=0; i<pop_.size(); i++){      //look at everything in the population vector
   c = pop_[i];                          //get the value of each chromosome 
   double fitness = c->get_fitness();   //get the fitness of each chromosome
   if (fitness > max){                    //find the maximum fitness value by checking every elements fitness value and replacing the current max if the next one is greated than the max
     max = fitness;
     cmax = c;
   }
  }
  return c;                             //return a copy of the chromosome with the highest fitness
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent()
{
  double sum = 0;               //set sum to 0
  for (int i=0; i<pop_.size(); i++){    //the fitness values for all of the chromosomes in the population vector
    double fitness = pop_[i]->get_fitness();    
    sum += fitness;             //add it to the sum 

  }
  double random_number = static_cast <double> (rand()) / static_cast <double> (RAND_MAX/sum);   //find a random number between 0 and the sum
  std::vector<Chromosome*> decending_chromosomes;   //make a new vector for the chromosomes in decending order
  for (int i=0; i<pop_.size(); i++){                //copy the pop_ vector to the decending chromosome vector
    decending_chromosomes.push_back(pop_[i]);
  }
  sort(decending_chromosomes.begin(), decending_chromosomes.end(), compare_fitness);    //sort the decending_chromosome vector in decending order
  double running_sum = 0;       //make a new sum
  for (int i=0; i<decending_chromosomes.size(); i++){   //do the roulette algorithm
    if(decending_chromosomes[i]->get_fitness()+running_sum > random_number){    //if the finess value of the chromosome is greater than that of the random number
      return decending_chromosomes[i];      //return the chromosome
    }
    running_sum += decending_chromosomes[i]->get_fitness();     //else add the next fitness value to the running_sum and compare that against the random number
  }
  return 0;
}
