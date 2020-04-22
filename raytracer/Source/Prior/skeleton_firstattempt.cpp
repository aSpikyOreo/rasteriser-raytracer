#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "limits"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;


#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

/* --------------------------------------- */
/* GLOBALS                                 */

vector<Triangle> triangles;

struct Intersection{
  vec4 position;
  float distance;
  int triangleIndex;
};

float m = std::numeric_limits<float>::max();

float focalLength = 1.0f;
vec4 cameraPos(SCREEN_WIDTH/2,SCREEN_HEIGHT/2,-3.0,1.0);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
bool CheckIntersectionInBounds(Triangle t, vec3 intersection);

bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection)
{
    bool intersectionExists = false;
    Intersection closest;
    for (uint i = 0; i < triangles.size(); i++){
      vec4 v0 = triangles[i].v0;
      vec4 v1 = triangles[i].v1;
      vec4 v2 = triangles[i].v2;

      vec3 e1 = vec3(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z); //v1-v0
      vec3 e2 = vec3(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z); //v2-v0
      vec3 b = vec3(start.x-v0.x, start.y-v0.y, start.z-v0.z); //start-v0

      vec3 d = vec3(dir.x,dir.y,dir.z);

      mat3 A(-d,e1,e2);
      vec3 x = glm::inverse(A)*b;

      if (CheckIntersectionInBounds(triangles[i], x)){
        //printf("intersection in bounds for triangle %d\n",i);
        if (!intersectionExists || closest.distance > x.x){
          intersectionExists = true;
          vec4 pos = vec4(x.x,x.y,x.z,1.0); //obviously needs correcting
          float dist = x.x;
          //closest = new Intersection();
          closest.position = pos;
          closest.distance = dist;
          closest.triangleIndex = i;
        }
      }
    }
    closestIntersection = closest;
    return intersectionExists;
}

bool CheckIntersectionInBounds(Triangle tri, vec3 intersection){
  float t = intersection.x;// - tri.v0.x;
  float u = intersection.y;// - tri.v0.y;
  float v = intersection.z;// - tri.v0.z;

  //if (t > 0) printf("t=%f; u=%f; v=%f\n",t,u,v );
  //return true;
  return (0 <= t && 0 < u && 0 < v && u+v < 1);
}


// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 colour;
// };

int main( int argc, char* argv[] )
{
  LoadTestModel(triangles);
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  while( NoQuitMessageSDL() )
    {
      Update();
      Draw(screen);
      SDL_Renderframe(screen);
      printf("done!\n");
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  int hatch = 0;
  for (int x = 0; x < screen->width; x++){
    for (int y = 0; y < screen->height; y++){
      //printf("x = %d, y = %d\n",x,y);
      vec4 start = cameraPos; //vec4(cameraPos.x+x-(screen->height/2),cameraPos.y+y-(screen->height/2),cameraPos.z,1.0);
      vec4 dir   = vec4(x-(screen->height/2),y-(screen->height/2),focalLength,1.0);
      //vec4 dir = vec4(0.0,0.0,1.0,1.0);
      Intersection intersection;
      if (ClosestIntersection(start,dir,triangles,intersection)){
        PutPixelSDL(screen,x,y,triangles[intersection.triangleIndex].color);
      } else {
        PutPixelSDL(screen,x,y,vec3(0.0,0.0,0.0));
      }
      hatch++;
      //if (hatch > 100) break;
    }
    //if (hatch > 100) break;
  }

  //vec3 colour(1.0,0.0,0.0);
  // for(int i=0; i<1000; i++)
  //   {
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  //   }
}

/*Place updates of parameters here*/
void Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  //std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}
