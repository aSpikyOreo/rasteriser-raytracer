#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 275
#define SCREEN_HEIGHT 275
#define FULLSCREEN_MODE false
#define USE_MATH_DEFINES


vector<Triangle> triangles;

float focalLength = 250.0f;
vec4 cameraPos(0.0f, 0.0f,-2.8f,1.0f);
vec4 cameraStart = cameraPos;

float a[16] = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
};

float yaw   = 0.0f;
float pitch = 0.0f;

struct Intersection
{
vec4 position;
float distance;
int triangleIndex;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void LoadTestModel( vector<Triangle>& triangles);
bool ClosestIntersection(vec4 start,vec4 dir,const vector<Triangle>& triangles,Intersection& closestIntersection );
float Magnitude(vec4 a, vec4 b);
mat3 RotateAllAxes(float Rx, float Ry, float Rz);
vec3 Intensity(vec3 power, float rad);

mat3 cameraRot = RotateAllAxes(0,0,0);
vec4 right_v  ( cameraRot[0][0], cameraRot[0][1], cameraRot[0][2], 1);
vec4 down_v   ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2], 1);
vec4 forward_v ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
vec4 lightPos( 0, -0.5, -0.7, 1.0 );
vec3 lightColor = 14.0f * vec3( 1, 1, 1 ); //defines power aka intensity of emitted light
vec3 DirectLight( const Intersection& i );
vec4 lightStart = lightPos;
vec3 indirectLight = 0.1f *vec3(1,1,1);




int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  printf("TEST\n");
  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }
cameraRot = RotateAllAxes(0.0f,yaw,0.0f);
  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}


vec3 CheckIntersection(Triangle t, vec4 s, vec4 d){
  //s = start
  //d = direction
  vec4 v0 = t.v0;
  vec4 v1 = t.v1;
  vec4 v2 = t.v2;

  vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
  vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
  vec3 b  = vec3(s.x-v0.x ,s.y-v0.y ,s.z-v0.z);

  //vec3 d_3 = vec3(d.x/d[3],d.y/d[3],d.z/d[3]);
  vec3 d_3 = vec3(d.x,d.y,d.z);


  mat3 A( -d_3, e1, e2 );
  vec3 x = glm::inverse( A ) * b;

  return x;
}

// vec3 colour(1.0,0.0,0.0);
// for(int i=0; i<1000; i++)
//   {
//     uint32_t x = rand() % screen->width;
//     uint32_t y = rand() % screen->height;
//     PutPixelSDL(screen, x, y, coloumat3 RotateAllAxes(float Rx, float Ry, float Rz)r);
//

float MAX = std::numeric_limits<float>::max();

mat3 RotateAllAxes(float Rx, float Ry, float Rz){
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

// vec3 cartesian_coordinate_rotate(vec3 coordinates, mat3 rotation_matrix){
//
// }
//
// mat4 RotateAroundY(vec3 position_to_rotate_at){
//   float a = position_to_rotate_at.x
//   float aaa[16] = {
//      1, 0, 0, 0,
//      5, 6, 7, 8,
//      9, 10, 11, 12,
//      13, 14, 15, 16
//   };
//   glm::mat4 bbb;
//   memcpy( glm::value_ptr( bbb ), aaa, sizeof( aaa ) );
//
// }
//
// vec4 right  ( cameraPos[0][0], cameraPos[0][1], cameraPos[0][2], 1);
// vec4 down   ( cameraPos[1][0], cameraPos[1][1], cameraPos[1][2], 1);
//

bool ClosestIntersection(vec4 start,vec4 dir,const vector<Triangle>& triangles,
Intersection& closestIntersection )
{
  closestIntersection.distance = MAX;
  bool intersected = false;
  //int TriangleNo = sizeof(triangles);
  for(int i = 0; i < triangles.size(); i++){
    vec3 x = CheckIntersection(triangles[i], start, dir);

    float t = x.x;
    float u = x.y;
    float v = x.z;

    if (0 <= u && 0 <= v && u + v <= 1 && 0 < t){
      //TODO; distance check {}
      if (t < closestIntersection.distance) {
        closestIntersection.position = start + t*dir;   //vec4(t,u,v,1);
        closestIntersection.distance = t;//Magnitude(closestIntersection.position,start);
        closestIntersection.triangleIndex = i;
        intersected = true;
      }
    }
  }
  return intersected;
}

//truncation error somewhere e.g. multiply float with int

vec3 Intensity(vec3 power, float rad, vec3 r, vec3 n){
  float area = (4 * M_PI * (rad*rad));
  vec3 B = power / area;
  return (B * max(glm::dot(glm::normalize(r),glm::normalize(n)), 0.0f));
  //absolute ???

  //x - lightPos | where lightPos is the light source
  //             | where x is the surface point

  // n           | where n is the normal pointing out of the surface
}


vec3 DirectLight(const Intersection &i){
  vec4 intersecting_point = i.position;
  vec3 r = glm::vec3(lightPos - intersecting_point);
  Intersection closestSurfaceToLight;

  float rad = glm::length(r);
  vec3 norm = glm::vec3(triangles[i.triangleIndex].normal);
  vec3 direct_light = Intensity(lightColor,rad,r,norm);


  bool closestPoint = ClosestIntersection(lightPos, intersecting_point-lightPos,triangles,closestSurfaceToLight);
  float distA = glm::length(i.position - lightPos);
  float distB = closestSurfaceToLight.distance;
  //printf("distA: %f vs distB: %f vs i_2.distance: %f\n",distA,distB,i_2.distance);
 if((distB < 0.999) && (distA > distB)){
   direct_light = vec3(0.0f,0.0f,0.0f);
 }
  return direct_light;
}



/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  LoadTestModel(triangles);
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 black = vec3(1.0,1.0,1.0);
  for (int x = 0; x < SCREEN_WIDTH; x++){
    for (int y = 0; y < SCREEN_HEIGHT; y++){
        //vec4 dir = vec4 (x - SCREEN_WIDTH/2,y - SCREEN_HEIGHT/2,focalLength,1);
        float local_x = x - SCREEN_WIDTH/2;
        float local_y = y - SCREEN_HEIGHT/2;
        float local_z = focalLength;

        vec4 global_x = local_x * right_v;
        vec4 global_y = local_y * down_v;
        vec4 global_z = local_z * forward_v;

        vec4 sum = global_x + global_y + global_z;
        vec4 dir = vec4 (sum.x, sum.y, sum.z, 1);
        //vec4 dir = vec4 (x - SCREEN_WIDTH/2,y - SCREEN_HEIGHT/2,focalLength,1);

        //dir *= ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
        Intersection i;
        bool intersected = ClosestIntersection(cameraPos,dir,triangles,i);
        if (intersected){
          //PutPixelSDL(screen,x,y,triangles[i.triangleIndex].color);
          PutPixelSDL(screen,x,y,triangles[i.triangleIndex].color*(DirectLight(i)+indirectLight));
        } else {
          PutPixelSDL(screen,x,y,black);
        }
    }
  }

  // vec3 colour(1.0,0.0,0.0);
  // for(int i=0; i<1000; i++)
  //   {
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  //   }
}

/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
        case SDLK_w:
            cameraPos += 0.01f * forward_v;
            break;
        case SDLK_s:
            cameraPos -= 0.01f * forward_v;
            break;
        case SDLK_a:
            cameraPos -= 0.01f * right_v;
            break;
        case SDLK_d:
            cameraPos += 0.01f * right_v;
            break;
        case SDLK_q:
            cameraPos -= 0.01f * down_v;
            break;
        case SDLK_e:
            cameraPos += 0.01f * down_v;
            break;
	      case SDLK_UP:    /* Move camera forward */
            pitch += 0.01;
            cameraRot = RotateAllAxes(pitch,0.0f,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2], 1 );
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2], 1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
            break;
	      case SDLK_DOWN:  /* Move camera backwards */
            pitch -= 0.01;
            cameraRot = RotateAllAxes(pitch,0.0f,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2], 1 );
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2], 1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
            break;
	      case SDLK_LEFT:	 /* Move camera left */
            yaw -= 0.01;
            cameraRot = RotateAllAxes(0.0f,yaw,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2], 1 );
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2], 1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
            break;
        case SDLK_RIGHT: /* Move camera right */
            yaw += 0.01;
            cameraRot = RotateAllAxes(0.0f,yaw,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2], 1 );
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2], 1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2], 1);
            //cameraPos.x += 0.01f;
            break;
        case SDLK_KP_8:
            lightPos -= 0.1f * down_v;
            break;
        case SDLK_KP_2:
            lightPos += 0.1f * down_v;
            break;
        case SDLK_KP_4:
            lightPos -= 0.1f * right_v;
            break;
        case SDLK_KP_6:
            lightPos += 0.1f * right_v;
            break;
	      case SDLK_ESCAPE:

      		/* Move camera quit */
      		return false;
        case SDLK_RETURN:
          lightPos = lightStart;
          cameraPos = cameraStart;
          break;
	      }
	  }
    }
  return true;
}
