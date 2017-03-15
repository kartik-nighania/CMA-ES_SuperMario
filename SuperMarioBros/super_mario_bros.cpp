/**
 * @file super_mario_bros.cpp
 * @author Marcus Edel
 *
 * Simple test programm to use the mlpack evolution methods.
 */

#include <mlpack/core.hpp>

#include <iostream>
#include <string>
#include <stdlib.h>

#include "parser.hpp"
#include "client.hpp"
#include "messages.hpp"


#include "link_gene.hpp"
#include "neuro_cmaes.hpp"
#include "neuron_gene.hpp"
#include "genome.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "random.hpp" 


/* #include <mlpack/methods/neuro_cmaes/neuro_cmaes.hpp>
#include <mlpack/methods/neuro_cmaes/random.hpp>
#include <mlpack/methods/neuro_cmaes/link_gene.hpp>
#include <mlpack/methods/neuro_cmaes/neuron_gene.hpp>
#include <mlpack/methods/neuro_cmaes/genome.hpp>
#include <mlpack/methods/neuro_cmaes/parameters.hpp>
#include <mlpack/methods/neuro_cmaes/utils.hpp>
*/




using namespace mlpack;
using namespace mlpack::neuro_cmaes;

class TaskSuperMarioBros
{
 public:
  /**
   * Create the super mario bros object.
   */
  TaskSuperMarioBros() { /* Nothing to do here */ }

  /**
   * Create the super mario bros object using the specified host and port.
   */
  TaskSuperMarioBros(const std::string& host,
                     const std::string& port) :
      host(host),
      port(port),
      success(false)
  {
     /* Nothing to do here */
  }

  /*
   * Send the next action to the connected server.
   *
   * client The client instance.
   */
  bool Action(const size_t action, client::Client& client)
  {
    try
    {
      if (action == 0)
      {
        client.Send(messages::JSONMessage(messages::PressRight()));
      }
      else if (action == 1)
      {
        client.Send(messages::JSONMessage(messages::PressLeft()));
      }
      else if (action == 2)
      {
        client.Send(messages::JSONMessage(messages::PressUp()));
      }
      else if (action == 3)
      {
        client.Send(messages::JSONMessage(messages::PressDown()));
      }
      else if (action == 4)
      {
        client.Send(messages::JSONMessage(messages::PressA()));
      }
    }
    catch (const std::exception& ex)
    {
      Log::Warn << ex.what() << std::endl;
      return false;
    }
    catch (...)
    {
      Log::Warn << "Send timeout." << std::endl;
      return false;
    }

    return true;
  }

  /*
   * Get the current game infromations from the connected server.
   *
   * @param client The client instance.
   */
  bool GameInfo(client::Client& client)
  {
    try
    {
      client.Send(messages::JSONMessage(messages::GameInfo()));

      std::string json;
      client.Receive(json);

      parser.Parse(json);
      parser.Tiles(tiles);

      parser.MarioPostion(marioPostionX, marioPostionY);
      parser.MarioLives(marioLives);
      parser.PlayerState(playerState);

    }
    catch (const std::exception& ex)
    {
      Log::Warn << ex.what() << std::endl;
      return false;
    }
    catch (...)
    {
      Log::Warn << "Receive timeout." << std::endl;
      return false;
    }

    return true;
  }

  /*
   * Fill the input vector with the screen infromations.
   *
   * @param input The vector used to store the game screen informations.
   */
  void DiscreteActuator(std::vector<double>& input)
  {
    for (size_t i = 0; i < tiles.n_elem; ++i)
    {
      input.push_back(tiles(i));
    }

    input.push_back(1.0);
  }

  /*
   * Connect to the server instance and reset the game state.
   *
   * @param client The client instance.
   */
  bool Reset(client::Client& client)
  {
    try
    {
      client::Client clientMaster;
      clientMaster.Connect(host, port);

      clientMaster.Send(messages::GetEndpoint());
      std::string json;
      clientMaster.Receive(json);

      parser.Parse(json);
      parser.Endpoint(hostEndpoint, portEndpoint);

      // Check if local balancer.
      if (hostEndpoint == "*") hostEndpoint = host;

      client.Connect(hostEndpoint, portEndpoint);
      client.Send(messages::JSONMessage(messages::ConfigSpeed("maximum")));
      client.Send(messages::JSONMessage(messages::ConfigDivisor(2)));
      client.Send(messages::JSONMessage(messages::PressRight()));
      client.Send(messages::JSONMessage(messages::GameReset()));
    }
    catch (const std::exception& ex)
    {
      Log::Warn << ex.what() << std::endl;
      return false;
    }
    catch (...)
    {
      Log::Warn << "Send timeout." << std::endl;
      return false;
    }

    return true;
  }

  // Whether task success or not.
  bool Success()
  {
    return success;
  }

  /*
   * Check if mario dies.
   */
  bool IsDead()
  {
    if (playerState == 11 || arma::accu(tiles) == 3)
    {
      return true;
    }

    return false;
  }

  /*
   * Evaluate the specified genome.
   *
   * @param genome Genome used for the evaluation process.
   */
  double EvalFitness(Genome& genome)
  {
    // Connect and reset game state.
    client::Client client;

    if(!Reset(client)) return 1;

    arma::mat tiles;
    size_t numSteps = 100000000;
    int maxMarioPositionX = 0;
    int sessionLives = 0;
    size_t stepCounter = 0;

    for (size_t step = 0; step < numSteps; ++step, ++stepCounter)
    {
      // Get the current game informations.
      if (!GameInfo(client)) continue;

      // Set the initial position and number of lives.
      if (step == 0)
      {
        maxMarioPositionX = marioPostionX;
        sessionLives = marioLives;
      }

      // Set network input.
      std::vector<double> input;
      DiscreteActuator(input);

      // Get network output.
      genome.Activate(input);
      std::vector<double> output;
      genome.Output(output);

      auto biggest_position = std::max_element(std::begin(output),
          std::end(output));
      size_t action = std::distance(std::begin(output), biggest_position);

      // Perform the action using the network output.
      if (!Action(action, client)) continue;

      // Check if mario dies.
      if (IsDead()) break;

      // Update marios position and reset the step counter.
      if (marioPostionX > maxMarioPositionX)
      {
        maxMarioPositionX = marioPostionX;
        stepCounter = 0;
      }

      // Abort if the marios x postion does not change in 80 steps.
      if (stepCounter >= 70) break;
    }

    // First level.
    if (maxMarioPositionX >= 3266)
    {
        success = true;
    }

    if (maxMarioPositionX > 0)
    {
      return 1 / double(maxMarioPositionX);
    }

    return 1;
  }

 private:
  //! Locally stored host name.
  std::string host;

  //! Locally stored port.
  std::string port;

  //! Locally stored endpoint host name.
  std::string hostEndpoint;

  //! Locally stored endpoint port.
  std::string portEndpoint;

  //! Locally stored parser instance.
  parser::Parser parser;

  //! Locally stored matrix that holds the tile information.
  arma::mat tiles;

  //! Locally stored x coordinate of mario.
  int marioPostionX;

  //! Locally stored y coordinate of mario.
  int marioPostionY;

  //! Locally stored number of lives.
  int marioLives;

  //! Locally stored player state.
  int playerState;

  //! Locally stored success indicator; set to true if task solved.
  bool success = false;
};


void setWeights(std::vector<LinkGene>& links, double const *x, int N)
{
     for(int i=0; i < N; i++) links[i].Weight(x[i]);
}

int main(int argc, char* argv[])
{
  mlpack::math::RandomSeed(1);

  std::string host(argv[1]);
  std::string port(argv[2]);

  TaskSuperMarioBros task(host, port);

  CMAES<double> evo;
  Parameters<double> params;

  params.stopMaxFunEvals= 50000;
  params.stopMaxIter= 10000;
  params.stStopFitness.flg = true;
  params.stStopFitness.val = 1/3266;
  params.logWarnings = true;
  params.lambda = 10;

   const int dim = 170*6 + 6*5;

  double xstart[dim];
  for(int i=0; i<dim; i++) xstart[i] = 0.5;

  double stddev[dim];
  for(int i=0; i<dim; i++) stddev[i] = 0.3;

   params.init(dim, xstart, stddev);

   double *arFunvals , *const*pop ;
   arFunvals = evo.init(params);

  // Set seed genome for the Super Mario Bros. task.
  ssize_t numInput = 170;
  ssize_t numOutput = 5;
  double fitness = -1;
  std::vector<NeuronGene> neuronGenes;
  std::vector<LinkGene> linkGenes;

  // Create number of input nodes.
  for (size_t i = 0; i < 169; ++i)
  {
    NeuronGene inputGene(i, INPUT, LINEAR, 0, 0, 0);
    neuronGenes.push_back(inputGene);
  }

  // Create bias node.
  NeuronGene biasGene(169, BIAS, LINEAR, 0, 0, 0);
  neuronGenes.push_back(biasGene);

  // Create output nodes.
  NeuronGene outputGene1(170, OUTPUT, SIGMOID, 1, 0, 0);
  neuronGenes.push_back(outputGene1);
  NeuronGene outputGene2(171, OUTPUT, SIGMOID, 1, 0, 0);
  neuronGenes.push_back(outputGene2);
  NeuronGene outputGene3(172, OUTPUT, SIGMOID, 1, 0, 0);
  neuronGenes.push_back(outputGene3);
  NeuronGene outputGene4(173, OUTPUT, SIGMOID, 1, 0, 0);
  neuronGenes.push_back(outputGene4);
  NeuronGene outputGene5(174, OUTPUT, SIGMOID, 1, 0, 0);
  neuronGenes.push_back(outputGene5);

  NeuronGene hiddenGene1(175, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene1);
  NeuronGene hiddenGene2(176, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene2);
  NeuronGene hiddenGene3(177, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene3);
  NeuronGene hiddenGene4(178, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene4);
  NeuronGene hiddenGene5(179, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene5);
  NeuronGene hiddenGene6(180, HIDDEN, SIGMOID, 0.5, 0, 0);
  neuronGenes.push_back(hiddenGene6);

  //connecting input to hidden nodes
  for(size_t j = 175 ; j <181 ; j++ )
  {
	  	for (size_t i = 0; i < 170; ++i)
	  {
	    LinkGene link(i,j,0);
	    linkGenes.push_back(link);
	  }
  }

  // connecting hidden to output nodes
   for(size_t j=175; j<181; j++)
   {
   	for (size_t i = 170; i < 175; ++i)
	  {
	     LinkGene link1(j, i, 0);
         linkGenes.push_back(link1);
	  }
   }



  // Instantiate initial seed genome.
  Genome neuralNet = Genome(neuronGenes, linkGenes, numInput, numOutput);
  neuralNet.SortLinkGenes();

 while(!evo.testForTermination() && !task.Success())
  {
    // Generate lambda new search points, sample population
    pop = evo.samplePopulation();

    // evaluate the new search points using fitness function from above
    for (int i = 0; i < evo.sampleSize(); ++i)
	{
		//wights and flush
		 neuralNet.Flush();

		 setWeights(linkGenes, pop[i], (int) evo.dimension());

		 arFunvals[i] = task.EvalFitness(neuralNet);

    }
    // update the search distribution used for sampleDistribution()
    evo.updateDistribution(arFunvals);
  }

  return 0;
}
