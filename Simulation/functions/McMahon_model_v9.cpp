// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;

// This is the model. It incorporates all the required variables to run replicates of the simulation. It can do so with a random start age for the individuals, and a random decay rate, if so desired.
//As currently written, this function will return stats to interrogate a single run if asked to do a single repeat, or if asked to do numerous repeat, a vector of the important stats required for these

// [[Rcpp::export]]
RcppExport SEXP
doSimulationAge (int individuals, float birthChance, float decayRate, int repeats, int startAge, int runFor, bool randomStartAge = false, bool randomDecayRate = false, bool timer = false)
{
  //There is also a timer included to allow for bench marking
  Timer benchmark;
  //If using the timer, there are a series of named benchmarking points included in the function outside the innner loop
  if (timer) benchmark.step ("Start");
  
  //We want to have sensible minima, as if either of these is zero, the model will not work
  int minimumStartAge = 10;
  float minimumDecayRate = 0.00001;
  
  //Conduct some checks
  if (randomStartAge && startAge < minimumStartAge)
  {
    Rcout << "startAge is smalller than the minimum randon startage - please modify code and retry";
    return 0;
  }
  
  if (randomDecayRate && decayRate < minimumDecayRate)
  {
    Rcout << "decayRate is smalller than the minimum randon startage - please modify code and retry";
    return 0;
  }
  
  // Data structure used to record the mean non zero value at equilibrium if we are doing multiple repeats - this is our "taphonomy value"
  // Use Rcpp vector here allowing us to return it directly
  Rcpp::NumericVector taphonomyVectorGlobal;
  
  // These data structures are here to record the outputs that can be used to investigate a single repeat.#
  // Use std vectors here for performance 
  std::vector<double> taphonomyVectorLocal;
  std::vector<double> differencesVector;
  std::vector<int> aliveVector;
  std::vector<int> decayingVector;
  std::vector<int> iterationVector;
  
  // These data structures exist to allow us to take a snapshot of the simulation at equilibrium
  std::vector<int> startAges;
  std::vector<int> states;
  std::vector<double> decayLevels;
  std::vector<double> decayRates;
  
  // Defines don't seem to work in Rcpp - use static variables instead
  static int dead = 0;
  static int alive = 1;
  static int decaying = 2;
  
  // Define my organism structure - the simulation comprises an organism list, which is of length *individuals* 
  struct organism
  {
    // Defaults in structure definitions fine from C++11
    float decayValue = 1.;
    float decayRate = 1.;
    int state = 0;
    int age = 0;
    int startAge = 0;
  };
  
  if (timer) benchmark.step ("Assign");
  
  //In this loop we run repeats of our simulation - either one or many
  for (int i = 0; i < repeats; i++)
  {
    //Internal data
    int iterations = 0;
    double meanNotDeadStored = 0.;
    double difference = 0.;
    bool warning = false;
    
    //Output some useful info to terminal
    Rcout << "\n\n***********Repeat " << i << "***********\n\n";
    Rcout << "\n\n Settings: \nIndividuals: " << individuals << "\nBirthChance: " << birthChance << "\nDecayRate: " << decayRate << "\nStartAge: " << startAge << "\nRunfor: " << runFor << "\nRandomStartAge: " << randomStartAge << "\nRandomDecayRate: " << randomDecayRate << "\n\n";
    
    // Make a vector of organisms - this is our list of living organisms, which can either be alive and of decreasing age, or yet to be spawned
    std::vector<organism> livingOrganisms(individuals);
    // Make a second list - this is our graveyard of dead & decaying organisms - they are added to this on death, and removed from it on complete decay
    std::vector<organism> graveyardOrganisms;
    
    // Can't use a variable local to this function to define this as default in the structure, therefore do it in a loop like this
    for (auto &o : livingOrganisms)
    {
      //If we are using a random start age, lets initialise this between the minimum, defined above, and the requested age which acts as the maximum
      if (randomStartAge) o.startAge = R::runif (minimumStartAge, startAge);
      else o.startAge = startAge;
      o.age = o.startAge;
      
      // We only need to set this here as it is not changed during the run - when reborn, and organism will keep its random decay rate
      // That said, does this lower the chances of a slowly decaying organism being reborn as there are fewer dead ones at any given time
      if (randomDecayRate) o.decayRate = R::runif (minimumDecayRate, decayRate);
      else o.decayRate = decayRate;
    }
    
    // Now we can do the actual simulation
    do
    {
      // Note here the need to pass by reference, or otherwise changes won't stick
      for (auto &o : livingOrganisms)
      {
        // If organism is not in existence, user set chance of spawning
        if (o.state == dead)
        {
          if (R::runif (0, 1) < birthChance)
          { 
            o.state++;
            if (randomDecayRate) o.decayRate = R::runif (minimumDecayRate, decayRate);
          }
        }
        // If it is is in existence, we have two choices, see below
        else if (o.state == alive)
        {
          //Either it has reached the end of its life. If so we move it to the graveyard list and reset it to not having spawned yet in this list
          if (o.age == 0)
          {
            //Kill the organism and place a copy in the graveyard list
            o.state++;
            organism deadOrganism = o;
            graveyardOrganisms.push_back (deadOrganism);
            // Reset the organism on the living list
            o.state = dead;
            o.decayValue = 1.;
            if (randomStartAge) o.startAge = R::runif (minimumStartAge, startAge);
            else o.startAge = startAge;
            o.age = o.startAge;
          }
          //Otherwise it remains alive and we just decrease its age
          else o.age--;
        }
      }
      
      //Now we deal with our graveyard -  first lets remove  completely decayed organisms
      //Given this is C++17, we can use a std::remove_if for an easy life
      //auto iterator = std::remove_if (graveyardOrganisms.begin (), graveyardOrganisms.end (), [] (organism o) { return o.decayValue <= 0; });
      //graveyardOrganisms.erase (iterator, graveyardOrganisms.end ());
      //Note that if we use C++20, this would be even simpler
      std::erase_if(graveyardOrganisms, [](organism o) { return o.decayValue <= 0; });
      
      // All remaining organisms are decaying, but not yet decayed - subtract the decay rate
      for (auto &o : graveyardOrganisms) o.decayValue = o.decayValue - o.decayRate;
      
      // Now let's work out the mean non zero value across the simulation
      int aliveCount = 0;
      int decayingCount = graveyardOrganisms.size ();
      double totalValue = 0.;
      
      //Lets work out the total decay value and alive plus decaying count
      for (auto o : livingOrganisms)
          if (o.state == alive)
            {
              aliveCount++;
              totalValue += o.decayValue;
            }
        for (auto o : graveyardOrganisms) totalValue += o.decayValue;
        
        int notDead = aliveCount + decayingCount;
        double meanNotDead = totalValue / static_cast<double> (notDead);
        //In some instances it is useful to track changes in this value through time
        difference = fabs (meanNotDead - meanNotDeadStored);
        meanNotDeadStored = meanNotDead;
        
        //Every ten iterations, lets record our stats
        if (iterations % 10 == 0)
        {
          taphonomyVectorLocal.push_back(meanNotDead);
          // Store data if one repeat
          if (repeats == 1)
          {
            differencesVector.push_back(difference);
            aliveVector.push_back(aliveCount);
            decayingVector.push_back(decayingCount);
            iterationVector.push_back(iterations);
          }
        }
        
        //It is useful to provide a warnign for saturation if list is small or birth chance is high- if alive list fills up, it will add artefacts to our results
        if (aliveCount > (individuals / 5) * 4 && warning == false)
        {
          Rcout << "Individuals list is >80% alive. Enlarge?";
          warning = true;
        }
        
        //Let's output updates every now and then
        if (iterations % 5000 == 0) Rcout << "Iteration " << iterations << "; alive count is " << aliveCount << "; decaying count is " << decayingCount << ".\n";
        iterations++;
    }
    while (iterations < runFor);
    
    if (timer) benchmark.step ("Simulation_complete");
    
    // With multiple repeats, we want to record a mean taphonomy value - sort out relevant data from vectors (sorry about all the benchmarking, there was originally a bottleneck here)
    std::vector<double>  taphonomyForMean (taphonomyVectorLocal.end () - 100, taphonomyVectorLocal.end ());
    if (timer) benchmark.step ("Analysis_01a");
    double total = std::reduce (taphonomyForMean.begin (), taphonomyForMean.end ());
    if (timer) benchmark.step ("Analysis_01b");
    double meanOfMeans = total / 100.;
    if (timer) benchmark.step ("Analysis_01c");
    //We will calculate a mean of the last 100 recorded means, so spanning (10*100) = 1000 iterations 
    taphonomyVectorGlobal.push_back (meanOfMeans);
    if (timer) benchmark.step ("Analysis_01d");
    
    //If we have a single repeat, this is set up to return a snapshot of the simulation for interrogation in R
    if (repeats == 1)
    {
      if (timer) benchmark.step ("Analysis_02a");
      for (auto &o : livingOrganisms)
      {
        startAges.push_back (o.startAge);
        states.push_back (o.state);
        decayLevels.push_back (o.decayValue);
        decayRates.push_back (o.decayRate);
      }
      for (auto &o : graveyardOrganisms)
      {
        startAges.push_back (o.startAge);
        states.push_back (o.state);
        decayLevels.push_back (o.decayValue);
        decayRates.push_back (o.decayRate);
      }
    }
    else taphonomyVectorLocal.clear(); 
    
    if (timer)   benchmark.step ("Analysis_02b");
  //Here ends the outer loop that controls the repeat number of the simulations
  }

  //We're all done with our repeats - now all  that is left is to return data
  if (repeats == 1)
  {
    if (timer)   benchmark.step ("Analysis_03a");
    //This can be quite slow, and we can return less data if we're in a standard run than one with either random start ages or decay rates - thus have a conditional return here
    if (!randomStartAge && !randomDecayRate)
    {
      Rcpp::DataFrame outDataFrame = DataFrame::create (Rcpp::Named ("differencesVector") = differencesVector, Rcpp::Named ("taphonomyVectorLocal") = taphonomyVectorLocal, Rcpp::Named ("aliveVector") = aliveVector, Rcpp::Named ("decayingVector") = decayingVector, Rcpp::Named ("iterationVector") = iterationVector);
      if (!timer) return outDataFrame;
    }
    else
    {
      Rcpp::DataFrame outDataFrame = DataFrame::create (Rcpp::Named ("startAges") = startAges, Rcpp::Named ("decayLevels") = decayLevels, Rcpp::Named ("states") = states, Rcpp::Named ("decayRates") = decayRates);
      if (!timer) return outDataFrame;
    }
    if (timer)   benchmark.step ("Analysis_03b");
  }
  //If multiple repeats, all we want to do is return our global mean taphonomy level across replicates list
  else if (!timer) return taphonomyVectorGlobal;
  
  //However, if we are benchmarking, we need to return our timings instead
  return Rcpp::NumericVector (benchmark);
}
//Fin.