#include <iostream>
#include "glm/glm.hpp"
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

// int read_glowTexture(SDL_Surface *surface){
//   SDL_Surface *glow_img = IMG_Load ("1027.jpg");
//
//   if( !glow_img){
//     printf("IMG_Load: %s\n", IMG_GetError() );
//     return 1;
//   }
//
//
//   //Draws image onto bolts
//   return 0;
//
//
// }



vec3 randomOffset(vec3 offsetLimit){
  float x_rand = -offsetLimit.x + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.x)));
  float y_rand = -offsetLimit.y + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.y)));
  float z_rand = -offsetLimit.z + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*offsetLimit.z)));
  vec3 result = vec3(x_rand,y_rand,z_rand);
  return result;
}

mat3 rotateSplit(float Rx, float Ry, float Rz){
  mat3 R = mat3 (1.0f);
  R[0][0] = cos(Ry) * cos(Rz);
  R[0][1] = cos(Ry) * sin(Rz);
  R[0][2] =  (-1) * sin(Ry);
  R[1][0] = (-1)*cos(Rx)*sin(Rz) + sin(Rx)*sin(Ry)*cos(Rz);
  R[1][1] =  cos(Rx)*cos(Rz) + sin(Rx)*sin(Ry)*sin(Rz);
  R[1][2] = sin(Rx)*cos(Ry);
  R[2][0] = sin(Rx)*sin(Rz) + cos(Rx)*sin(Ry)*cos(Rz);
  R[2][1] = (-1)*sin(Rx)*cos(Rz) + cos(Rx)*sin(Ry)*sin(Rz);
  R[2][2] = cos(Rx)*cos(Ry);

  return R;
}

//algorithm that produces info required to draw the lightning

Lightning validateSegment(Lightning segment_to_check){
  vec3 start = segment_to_check.start;
  vec3 end = segment_to_check.end;

  if(start.x > 1) segment_to_check.start.x = 1;
  if(start.y > 1) segment_to_check.start.y = 1;
  if(start.z > 1) segment_to_check.start.z = 1;
  if(start.x < (-1)) segment_to_check.start.x = (-1);
  if(start.x < (-1)) segment_to_check.start.y = (-1);
  if(start.x < (-1)) segment_to_check.start.z = (-1);
  if(end.x > 1) segment_to_check.end.x = 1;
  if(end.y > 1) segment_to_check.end.y = 1;
  if(end.z > 1) segment_to_check.end.z = 1;
  if(end.x < (-1)) segment_to_check.end.x = (-1);
  if(end.x < (-1)) segment_to_check.end.y = (-1);
  if(end.x < (-1)) segment_to_check.end.z = (-1);

  return segment_to_check;
}

//sets a random integer between 1 and N
int get_random_index(int N){
  int index;
  srand(time(0));
  index = int(rand() % N) + 1;
  return index;
}

//sets a random co-ordinate to be an end-point for the lightning's path
vec3 setGoal(){
  vec3 limit = 0.99f*vec3(1,1,1);
  vec3 possibleGoal = randomOffset(limit);
  float val = 1.0f;
// 1 is X
// 2 is Y
// 3 is Z
int chosenPlane = get_random_index(3);
// 1 is -ve
// 2 is +ve
int parity = get_random_index(2);

if(parity == 1) val = -1.0f;


switch (chosenPlane) {
  case 1: possibleGoal.x = val;
  case 2: possibleGoal.y = val;
  case 3: possibleGoal.z = val;
}

  return possibleGoal;
}


vector<Lightning> drawOrigin(vec3 origin){
  float theta, rho, phi,phase;
  vector<Lightning> endSet;
  vec3 originOffset = 1.11f*vec3(1-phase,1+phase,1-phase);
  for(int l = 0; l < 10; l++){
    vec3 end = cos(rho)*rotateSplit(theta, rho, phi) * randomOffset(2.0f*originOffset*sin(theta));
    Lightning rotatingBolt = Lightning(originOffset,end);
    float lengthScale = 0.7f;
    vec3 currentMidpoint, direction, splitTransform;

    currentMidpoint = 0.5f*(end +origin);
    vec3 p1 = glm::normalize(end - origin);
    vec3 currentOffset = originOffset;
    //vec3 p2 = randomOffset(vec3(1,1,1));

    currentMidpoint += glm::cross(p1,vec3(0,0.98,0)) *randomOffset(currentOffset);
    theta+=0.01;
    rho -= 0.01;
    phi += 0.01;
    phase +=0.03;
    endSet.push_back(Lightning(originOffset,end));
  }
  return endSet;
}



vector<Lightning> ProduceLightning( vector<Lightning> currentBoltSegments, vec3 norm){
  vec3 point_of_origin = vec3(0.01,-0.62,-0.6);
  vec3 endGoal = setGoal();
  Lightning firstBolt = Lightning(point_of_origin, endGoal);
  currentBoltSegments.push_back(firstBolt);
  float lengthScale = 0.7f;
  vec3 maxOffset = 0.7f*vec3(1,1,1);
  vec3 currentMidpoint, direction, splitTransform;
  currentMidpoint = 0.5f*(endGoal +point_of_origin);
  vec3 p1 = glm::normalize(endGoal - point_of_origin);
  vec3 currentOffset = maxOffset;
  //vec3 p2 = randomOffset(vec3(1,1,1));

  currentMidpoint += glm::cross(p1,norm) *randomOffset(currentOffset);
  int N = 50;
  // for(int g = 0; g < 2; g++){
    // for(uint32_t i = 0; i < N; i++){
    //   currentBolt = currentBoltSegments[i];
    //   currentBoltSegments.erase(currentBoltSegments.begin()+i);
    // }

    endGoal = setGoal();
    Lightning currentBolt = Lightning(currentMidpoint,endGoal);
    currentMidpoint = 0.5f*(currentBolt.end + currentBolt.start);
    p1 = glm::normalize(currentBolt.end - currentBolt.start);
    //vec3 p2 = randomOffset(vec3(1,1,1));

    currentMidpoint += glm::cross(p1,norm) *randomOffset(currentOffset);
    ////////////////////////////////////////////////////////////////////////////
    direction = currentMidpoint - currentBolt.start;
    splitTransform = lengthScale*rotateSplit(0.13,0.39,0.27)*direction;
    /////////////////////////////////////////////////////////////////////////////
    Lightning split_segment_one = Lightning(currentBolt.start, currentMidpoint);
    Lightning split_segment_two = Lightning(currentMidpoint, currentBolt.end);
    Lightning split_segment_three = Lightning(currentMidpoint,splitTransform);
    split_segment_one = validateSegment(split_segment_one);
    split_segment_two = validateSegment(split_segment_two);
    split_segment_three =validateSegment(split_segment_three);
    currentBoltSegments.push_back(split_segment_one);
    currentBoltSegments.push_back(split_segment_two);
    currentBoltSegments.push_back(split_segment_three);
  //  }
  //
   currentOffset /= 2;
   //currentBolt = split_segment_three;
  //
  //   }

    return currentBoltSegments;
  }
