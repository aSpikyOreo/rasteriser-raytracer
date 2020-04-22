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


class Lightning{
  public:
    vec3 start;
    vec3 end;


    Lightning(glm::vec3 start, glm::vec3 end)
     : start(start), end(end)
     {

     }
};

vec3 randomOffset(vec3 offsetLimit){
  float x_rand = -offsetLimit.x + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.x)));
  float y_rand = -offsetLimit.y + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.y)));
  float z_rand = -offsetLimit.z + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.z)));
  vec3 result = vec3(x_rand,y_rand,z_rand);
  return result;
}

//algorithm that produces info required to draw the lightning
void ProduceLightning( vector<Lightning> currentBoltSegments, vec3 norm){
  Lightning firstBolt = Lightning(vec3(0,0,0), vec3(1,-0.4,0.6));
  currentBoltSegments.push_back(firstBolt);
  vec3 maxOffset = 0.7f*vec3(1,1,1);
  vec3 currentMidpoint;
  Lightning currentBolt = firstBolt;
  vec3 currentOffset = maxOffset;
  for(int g = 0; g < 5; g++){
    for(uint32_t i = 0; i < currentBoltSegments.size(); i++){
      currentBolt = currentBoltSegments[i];
      currentBoltSegments.erase(currentBoltSegments.begin()+i);

      currentMidpoint = 0.5f*(currentBolt.end + currentBolt.start);
      vec3 p1 = glm::normalize(currentBolt.end - currentBolt.start);
      //vec3 p2 = randomOffset(vec3(1,1,1));

      currentMidpoint += glm::cross(p1,norm) *randomOffset(currentOffset);
      Lightning split_segment_one = Lightning(currentBolt.start, currentMidpoint);
      currentBoltSegments.insert(currentBoltSegments.begin()+i, split_segment_one);
      Lightning split_segment_two = Lightning(currentMidpoint, currentBolt.end);
      currentBoltSegments.insert(currentBoltSegments.begin()+(i+1), split_segment_two);
    }
  currentOffset /= 2;

    }
    for(uint32_t b = 0; b < currentBoltSegments.size(); b++){
    printf("Bolt start: (%f ,%f ,%f) --> Bolt-end: (%f, %f, %f)",currentBoltSegments[b].start.x,
      currentBoltSegments[b].start.y,currentBoltSegments[b].start.z, currentBoltSegments[b].end.x,
        currentBoltSegments[b].end.y,currentBoltSegments[b].end.z);
    }
  }



int main( int argc, char* argv[] )
{
  vector<Lightning> currentBoltSegments;

  while ( 1 )
    {

      ProduceLightning(currentBoltSegments, vec3(1,1,1));

    }
  return 0;
}
