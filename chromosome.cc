/*
 * Implementation for Chromosome class
 */

#include <algorithm>
#include <cassert>
#include <random>
#include <iostream>

#include "chromosome.hh"
#include "cities.hh"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Generate a completely random permutation from a list of cities
Chromosome::Chromosome(const Cities* cities_ptr)
  : cities_ptr_(cities_ptr),
    order_(random_permutation(cities_ptr->size())),
    generator_(rand())
{
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Clean up as necessary
Chromosome::~Chromosome()
{
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Perform a single mutation on this chromosome
void
Chromosome::mutate()
{
  assert(is_valid());
  int firstIndex = generator_() % order_.size();
  int secondIndex;
  do {
    secondIndex = generator_() % order_.size();
  } while (secondIndex == firstIndex && order_.size() != 1);
  int tempElement = order_[firstIndex];
  order_[firstIndex] = order_[secondIndex];
  order_[secondIndex] = tempElement;
}

//////////////////////////////////////////////////////////////////////////////
// Return a pair of offsprings by recombining with another chromosome
// Note: this method allocates memory for the new offsprings
std::pair<Chromosome*, Chromosome*>
Chromosome::recombine(const Chromosome* other)
{
  assert(is_valid());
  assert(other->is_valid());

  // simple chromosomes are not complex enough to recombine
  if(order_.size() < 3) {
    return make_pair(this->clone(), other->clone());
  }

  // picks an index near the middle to split for crossover
  unsigned secondIndex = order_.size() / 2;
  
  pair<Chromosome*, Chromosome*> returnPair(create_crossover_child(this, other, 0, secondIndex),
					    create_crossover_child(other, this, 0, secondIndex));
  return returnPair;			    
}

//////////////////////////////////////////////////////////////////////////////
// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e),
// and all the other values in the same order as in p2.
Chromosome*
Chromosome::create_crossover_child(const Chromosome* p1, const Chromosome* p2,
                                   unsigned b, unsigned e) const
{
  Chromosome* child = p1->clone();

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b and i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
      }
      assert(j < p2->order_.size());
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
}

// Return a positive fitness value, with higher numbers representing
// fitter solutions (shorter total-city traversal path).
// uses the inverse function on distance to obtain fitness
double
Chromosome::get_fitness() const
{
  // avoid dividing by zero; return value technically arbitrary but
  // having two cities on top of each other isn't an expected edge case
  if (calculate_total_distance() == 0) {
    return 1e10;
  } else {
    return 1. / calculate_total_distance();
  }
}

// A chromsome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
// We implement this check with a sort, which is a bit inefficient, but simple
bool
Chromosome::is_valid() const
{
  // standard perm in the form {0, 1, 2, ...}
  Cities::permutation_t standardPermutation;
  standardPermutation.reserve(order_.size());
  standardPermutation.resize(order_.size());
  iota(standardPermutation.begin(), standardPermutation.end(), 0);

  // order sorted ascending
  Cities::permutation_t sortedOrder(order_);
  sort(sortedOrder.begin(), sortedOrder.end());

  // checks order has correct number of elements and doesn't skip
  // whole numbers.
  if (sortedOrder == standardPermutation) {
    return true;
  } else {
    return false;
  }			  
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.
bool
Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const
{
  if (find(order_.cbegin() + begin, order_.cbegin() + end, value) == order_.cbegin() + end) {
    return false;
  } else {
    return true;
  }
}
