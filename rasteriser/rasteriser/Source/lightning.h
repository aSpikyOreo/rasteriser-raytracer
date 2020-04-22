#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <vector>
#include<ctime>
#include<cstdlib>

using namespace std;
using glm::vec3;
using glm::vec2;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;

// describes a particular bolt of lightning




class Lightning{
  public:
    vec3 start;
    vec3 end;


    Lightning(glm::vec3 start, glm::vec3 end)
     : start(start), end(end)
     {

     }
};


vector<Lightning> boltSegments;

void randomOffset(vec3 offsetLimit){
  float x_rand = -offsetLimit.x + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.x)));
  float y_rand = -offsetLimit.y + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.y)));
  float z_rand = -offsetLimit.z + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.z)));
  vec3 result = vec3(x_rand,y_rand,z_rand);
  return result;
}


void ProduceLightning(std::vector<Lightning> currentBoltSegments, int n){
  currentBoltSegments.push_back(new Lightning(vec3(0,0,0), vec3(1,-0.4,0.6)));
  vec3 maxOffset = 0.7f*vec3(1,1,1);
  vec3 currentMidpoint;
  vec3 currentNormal;
  Lightning currentBolt;
  vec3 currentOffset = maxOffset;
  for(int g = 0; g < 5; g++){
    for(int i = 0; i < currentBoltSegments.size(); i++){
      currentBolt = currentBoltSegments[i];
      currentBoltSegments.erase(i);
      currentMidpoint = 0.5f*(currentBolt.end + currentBolt.start);
      vec3 p1 = glm::normalize(currentBolt.end - currentBolt.start);
      vec3 p2 = randomOffset(vec3(1,1,1));
      currentMidpoint += glm::perp(p1,p2) *randomOffset(currentOffset);

      currentBoltSegments.insert(i, new Lightning(currentBolt.start, currentMidpoint));
      currentBoltSegments.insert(i+1, new Lightning(currentMidpoint, currentBolt.end));
    }
  currentOffset /= 2;

    }
  }
