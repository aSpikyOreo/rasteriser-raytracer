#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include "Lightning.h"

using namespace std;
using glm::vec3;
using glm::vec2;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;

SDL_Event event;

#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define x_max SCREEN_WIDTH
#define y_max SCREEN_HEIGHT
#define z_max 150
#define FULLSCREEN_MODE false

#define MAX_BOLTS_ONSCREEN 90

vector<Triangle> triangles;
vector<Triangle> clippedTriangles;
vector<Lightning> boltSegments;
vector<Lightning> mirrorBolts; //attempt at adding the illusion of a 3rd dimension to the line
vector<Lightning> boltAvg;
vector<Lightning> originBolt;

// outlines whether the maximum bolts allowed on-screen
// has been exceeded or not
bool boltSwitch = true;

float focalLength = 250.0f;
vec4 cameraPos(0.0f,0.0f,-2.8f,1.0f);
vec4 lightPos(0.0, -0.45, -0.55,1.0f);
vec4 lightStart = lightPos;
vec3 lightPower = 14.0f*vec3( 1, 1, 1 );
vec3 indirectIntensity = 0.5f*vec3( 1, 1, 1 );
vec4 cameraStart = cameraPos;
float pitch = 0.0f;
float yaw = 0.0f;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float textureBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];


struct Pixel
{
int x;
int y;
float zinv;
vec3 illumination;
vec4 pos3d;
int u,v; //texture points
};

struct Vertex{
  vec4 pos;
  vec4 originalPos;
  vec4 pos_to_clip;
  int triangleIndex;
  int vertexIndex;
};

struct Texture{
  string textureName;
  vec3 normal;

  Vertex vertexA;
  Vertex vertexB;
  Vertex vertexC;
  vec3 color;
};


Texture boltTexture;
Texture wall_Texture;

typedef int OutCode;


const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

vec3 Intensity(vec3 power, float rad, vec3 r, vec3 n);
vec3 DirectLight(const Pixel &p);
void VertexShader(const Vertex& v, Pixel& p );
void PixelShader( const Pixel& p );
mat4 RotateAllAxes(float Rx, float Ry, float Rz);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void DrawLineSDL(screen* surface, Pixel a, Pixel b, vec3 color);
void DrawPolygon( const vector<Vertex>& vertices, screen* screen );
Vertex viewFrustrum(Pixel p,Vertex v);
void clipPolygon(Vertex v, vec4 P);
void clipAlgorithm(Pixel p, vec4 P);
void PixelsOnLine(Pixel a, Pixel b, vector<Pixel>& line);
void ComputePolygonRows(const vector<Pixel>& vertexPixels,
vector<Pixel>& leftPixels,vector<Pixel>& rightPixels );
Vertex validateVertex(Vertex v, vec2 x_limit,vec2 y_limit,vec2 z_limit);
void DrawPolygonRows(const vector<Pixel>& leftPixels,
const vector<Pixel>& rightPixels, screen* s, vec3 col);
vector<Lightning> allocateBolts(vector<Lightning> current_bolt_streaks);
Texture fillTextureBuffer(char* img_path);

/*                               ------------------------                       */
OutCode CalculateOutCode(float x, float y);
//essentially the easiest of the forms of clipping available
void Cohen_Suther_Clipping(screen* surface, Pixel a, Pixel b);


/*-------------------------------------------------------------------------------*/
mat4 cameraRot = RotateAllAxes(0,0,0);
vec4 currentNormal;
vec3 currentColour;
vec3 currentReflectance;


vec4 right_v  ( 1, 0, 0,1);
vec4 down_v   ( 0, 1, 0,1);
vec4 forward_v ( 0, 0, 1,1);



//
// Texture fillTextureBuffer(){
//   Texture newTexture;
//   if (SDL_LoadBMP("side_walls.bmp") == NULL) {
//     // Unrecoverable error, exit here.
//     printf("SDL_LoadBMP failed: %s\n", SDL_GetError());
// }
//
//   SDL_Surface* img = nullptr;
//   SDL_Texture *txt = nullptr;
//   img = SDL_LoadBMP("side_walls.bmp");
//   if(img != nullptr){
//     txt = SDL_CreateTexture(renderer,img);
//
//   }
//   uint32_t img_width = img->w;
//   uint32_t img_height = img->h;
//   vec3 img_col;
//   for(int jj = 0; jj < SCREEN_HEIGHT; jj++){
//     for(int ii = 0; ii < SCREEN_WIDTH; ii++){
//       img_col = vec3(img->format->Rmask,img->format->Gmask,img->format->Bmask);
//     }
//   }
//   newTexture.color = img_col;
//   return newTexture;
//
// }



OutCode CalculateOutCode(float x, float y){
  OutCode code;

  code = INSIDE;

  if (x < 0){        // to the left of clip window
		code |= LEFT;
  }
	else if (x > SCREEN_WIDTH){      // to the right of clip window
		code |= RIGHT;
  }
	if (y < 0){           // below the clip window
		code |= BOTTOM;
  }
	else if (y > SCREEN_HEIGHT){      // above the clip window
		code |= TOP;
  }

	return code;
}

void Cohen_Suther_Clipping(screen* surface,Pixel a, Pixel b){
  int x_min = 50;
  int y_min = 50;
  // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
	OutCode outcode0 = CalculateOutCode(a.x, a.y);
	OutCode outcode1 = CalculateOutCode(b.x, b.y);
	bool accept = false;

	while (true) {
		if (!(outcode0 | outcode1)) {
			// bitwise OR is 0: both points inside window; trivially accept and exit loop
			accept = true;
			break;
		} else if (outcode0 & outcode1) {
			// bitwise AND is not 0: both points share an outside zone (LEFT, RIGHT, TOP,
			// or BOTTOM), so both must be outside window; exit loop (accept is false)
			break;
		} else {
			// failed both tests, so calculate the line segment to clip
			// from an outside point to an intersection with clip edge
			double x, y;

			// At least one endpoint is outside the clip rectangle; pick it.
			OutCode outcodeOut = outcode0 ? outcode0 : outcode1;

			if (outcodeOut & TOP) {           // point is above the clip window
				x = a.x + (b.x - a.x) * (y_max - a.y) / (b.y - a.y);
				y = y_max;
			} else if (outcodeOut & BOTTOM) { // point is below the clip window
				x = a.x + (b.x - a.x) * (y_min *a.y) / (b.y - a.y);
				y = y_min;
			} else if (outcodeOut & RIGHT) {  // point is to the right of clip window
				y = a.y + (b.y - a.y) * (x_max - a.x) / (b.x - a.x);
				x = x_max;
			} else if (outcodeOut & LEFT) {   // point is to the left of clip window
				y = a.y + (b.y - a.y) * (x_min - a.x) / (b.x - a.x);
				x = x_min;
			}

			// Now we move outside point to intersection point to clip
			// and get ready for next pass.
			if (outcodeOut == outcode0) {
				a.x = x;
				a.y = y;
				outcode0 = CalculateOutCode(a.x, a.y);
			} else {
				b.x = x;
				b.y = y;
				outcode1 = CalculateOutCode(b.x, b.y);
			}
		}
	}
	if (accept) {
		// Following functions are left for implementation by user based on
		// their platform (OpenGL/graphics.h etc.)

    //TODO:
		DrawLineSDL(surface, a,b,currentColour);
	}
}



mat4 PerspectiveTransform(float Rx, float Ry, float Rz){
  mat4 R = mat4 (1.0f);
  R[0][0] = cos(Ry) * cos(Rz);
  R[1][0] = cos(Ry) * sin(Rz);
  R[2][0] =  (-1) * sin(Ry);
  R[0][3] = 0;
  R[0][1] = (-1)*cos(Rx)*sin(Rz) + sin(Rx)*sin(Ry)*cos(Rz);
  R[1][1] =  cos(Rx)*cos(Rz) + sin(Rx)*sin(Ry)*sin(Rz);
  R[2][1] = sin(Rx)*cos(Ry);
  R[1][3] = 0;
  R[0][2] = sin(Rx)*sin(Rz) + cos(Rx)*sin(Ry)*cos(Rz);
  R[1][2] = (-1)*sin(Rx)*cos(Rz) + cos(Rx)*sin(Ry)*sin(Rz);
  R[2][2] = cos(Rx)*cos(Ry);
  R[2][3] =  (1/focalLength);
  R[3][0] = (-1)*cameraPos.x;
  R[3][1] = (-1)*cameraPos.y;
  R[3][2]= (-1)*cameraPos.z;
  R[3][3] = 0;

  return R;
}

mat4 ViewportTransform(){
  mat4 R = mat4 (1.0f);
  R[0][0] = 0.5;
  R[1][0] = 0;
  R[2][0] =  0;
  R[0][3] = 0;
  R[0][1] = 0;
  R[1][1] =  0.5;
  R[2][1] = 0;
  R[1][3] = 0;
  R[0][2] = 0;
  R[1][2] = 0;
  R[2][2] = 0.5;
  R[2][3] =  0;
  R[3][0] = 0.5;
  R[3][1] = 0.5;
  R[3][2]= 0.5;
  R[3][3] = 1;

  return R;
}




mat4 RotateAllAxes(float Rx, float Ry, float Rz){
  mat4 R = mat4 (1.0f);
  R[0][0] = cos(Ry) * cos(Rz);
  R[1][0] = cos(Ry) * sin(Rz);
  R[2][0] =  (-1) * sin(Ry);
  R[0][3] = 0;
  R[0][1] = (-1)*cos(Rx)*sin(Rz) + sin(Rx)*sin(Ry)*cos(Rz);
  R[1][1] =  cos(Rx)*cos(Rz) + sin(Rx)*sin(Ry)*sin(Rz);
  R[2][1] = sin(Rx)*cos(Ry);
  R[1][3] = 0;
  R[0][2] = sin(Rx)*sin(Rz) + cos(Rx)*sin(Ry)*cos(Rz);
  R[1][2] = (-1)*sin(Rx)*cos(Rz) + cos(Rx)*sin(Ry)*sin(Rz);
  R[2][2] = cos(Rx)*cos(Ry);
  R[2][3] = 0;
  R[3][0] = (-1)*cameraPos.x;
  R[3][1] = (-1)*cameraPos.y;
  R[3][2]= (-1)*cameraPos.z;
  R[3][3] = 1;

  return R;
}



















void PixelShader(screen* surface, const Pixel& p, vec3 color ){
  int x = p.x;
  int y = p.y;
  //Texture wall_tex = fillTextureBuffer("side_walls.bmp");
  if(p.zinv > depthBuffer[y][x]){
    depthBuffer[y][x] = p.zinv;
     PutPixelSDL(surface,x, y, currentColour*(DirectLight(p) + indirectIntensity));
  }
  boltSegments = allocateBolts(boltSegments);
  //originBolt = allocateBolts(originBolt);
}



//originalPos has camera co-coordinates
//pos has homogeneous coordinates
void VertexShader(const Vertex& v, Pixel& p,int t ){
   Vertex v_new = v;
   vec4 pos = v.pos;
   vec4 some_pos;

   p.zinv = 1 / pos.z;
   //projection onto screen space
   p.x = int(focalLength *pos.x*p.zinv + SCREEN_WIDTH/2);
   p.y = int(focalLength *pos.y*p.zinv + SCREEN_HEIGHT/2);
   some_pos = pos * PerspectiveTransform(pitch,yaw,0.0f);
   v_new.pos_to_clip = some_pos;
   v_new = viewFrustrum(p,v_new);

   p.pos3d = v.originalPos * p.zinv;

    //could possibly be originalPos rather than pos
    // if it doesn't work, make this alteration

}

Vertex viewFrustrum(Pixel p, Vertex v){
  vec4 _inClipSpace_ = v.pos_to_clip;
  float z_far = v.pos_to_clip.z;
  //normalising all co-ordinates to produce our NDC space
  _inClipSpace_.x /= _inClipSpace_.w;
  _inClipSpace_.y /= _inClipSpace_.w;
  _inClipSpace_.z /= _inClipSpace_.w;

  vec3 cubeNDC = vec3(_inClipSpace_.x, _inClipSpace_.y,_inClipSpace_.z);

   vec2 clipped_Xplane = vec2(_inClipSpace_.w * x_max * (-1), _inClipSpace_.w * x_max);
   vec2 clipped_Yplane = vec2(_inClipSpace_.w * y_max * (-1), _inClipSpace_.w * y_max);
   vec2 clipped_Zplane = vec2(_inClipSpace_.w * z_max * (-1), _inClipSpace_.w * z_max);

   // a = p.triangleIndex
   // v (ref to triangle it resides on) <---attached---> new Triangle[a + 1]
   //


   // p.in_clipping_space = false;
   // if(cubeNDC.x > clipped_Xplane.x && cubeNDC.x < clipped_Xplane.y){
   //   if(cubeNDC.y > clipped_Yplane.x && cubeNDC.y < clipped_Yplane.y){
   //     p.in_clipping_space = true;
   //   }
   // }
   //
   // if(cubeNDC)


//   //update to screen_space
//   _inClipSpace_.x = (_inClipSpace_.x + 1) * 0.5f * (x_max - 1);
//   _inClipSpace_.y = (_inClipSpace_.y + 1) * 0.5f * (y_max - 1);
// //  vec2 clippingLimits = vec2(_inClipSpace_.x, _inClipSpace_.y);

  //map clipped co-ordinates back
  // [ u , v , f , 1]

  vec4 clippedCoordinates = ViewportTransform() * vec4(_inClipSpace_.x,
                            _inClipSpace_.y, _inClipSpace_.z,1);

  float z_near = _inClipSpace_.z;
  float new_zinv = (z_far - z_near)*clippedCoordinates.z + (z_far-z_near)*1
          + z_near*clippedCoordinates.z +z_near;


  v.originalPos.x = clippedCoordinates.x * SCREEN_WIDTH;
  v.originalPos.y = clippedCoordinates.y * SCREEN_HEIGHT;
  p.x = v.originalPos.x;
  p.y = v.originalPos.y;

  v = validateVertex(v, clipped_Xplane, clipped_Yplane, clipped_Zplane);
  return v;
}

  /*
  Clipping takes place
  */
  //clipPolygon(v, clippedCoordinates);



  // void clipPolygon(Vertex v, vec4 P){
  //   float d1,d2,t;
  //   vec4 n, intersect_point;
  //   vec4 q1 = triangles[v.triangleIndex].v0;
  //   vec4 q2 = triangles[v.triangleIndex].v1;
  //   //float on_plane_clause = 0.00001f;
  //   n = currentNormal;
  //   d1 = glm::dot((q1 - P ), n);
  //   d2 = glm::dot((q2 - P ), n);
  //
  //   t = d1 / (d1 - d2);
  //   intersect_point = q1 + t*(q2 - q1);
  //   //Line completely inside
  //   if( ( (d1 >= 0) && (d2 >= 0)) || ((d1 >= 0) && (d2 >= 0))) {
  //   }
  //
  //   //Line completely outside
  //   if( ((d1 <= 0) && (d2 < 0)) || ((d1 < 0) && (d2 <= 0))){
  //
  //
  //   }
  //
  //   //Q1 inside and Q2 outside
  //   if((d1 > 0) && (d2 < 0)){
  //     intersect_point = q1 + t*(q2 - q1);
  //     q2 = intersect_point;
  //     //Q1_I lies on the in side
  //     //I_Q2 lies on the out side
  //   }
  //
  //   //Q2 inside and Q1 outside
  //   if((d1 < 0) && (d2 > 0)){
  //     intersect_point = q1 + t*(q2 - q1);
  //     q1 = intersect_point;
  //     //Q1_I lies on the out side
  //     //I_Q2 lies on the in side
  //   }
  //
  //   triangles[v.triangleIndex].v0 = q1;
  //   triangles[v.triangleIndex].v1 = q2;
  //
  // }





Vertex validateVertex(Vertex v, vec2 x_limit,vec2 y_limit,vec2 z_limit){
  vec4 clippedPos = v.originalPos;
  if(clippedPos.x < x_limit.x) clippedPos.x = x_limit.x;
  if(clippedPos.x > x_limit.y) clippedPos.x = x_limit.y;
  if(clippedPos.y < y_limit.x) clippedPos.y = y_limit.x;
  if(clippedPos.y > y_limit.y) clippedPos.y = y_limit.y;
  if(clippedPos.z < z_limit.x) clippedPos.z = z_limit.x;
  if(clippedPos.z > z_limit.y) clippedPos.z = z_limit.y;

  if(v.vertexIndex == 0) triangles[v.triangleIndex].v0 = clippedPos;
  if(v.vertexIndex == 1) triangles[v.triangleIndex].v1 = clippedPos;
  if(v.vertexIndex == 2) triangles[v.triangleIndex].v2 = clippedPos;


  v.originalPos = clippedPos;

  return v;
}


vec3 Intensity(vec3 power, float rad, vec3 r, vec3 n){
  float area = (4 * M_PI * (rad*rad));
  vec3 B = power / area;
  return (B * max(glm::dot(glm::normalize(n),glm::normalize(r)), 0.0f));
}


void DrawLineSDL(screen* surface, Pixel a, Pixel b, vec3 color){
  vector<Pixel> line(1);
  PixelsOnLine(a,b,line);
  int x,y;
  for (size_t i = 0; i < line.size(); i++){
    x = line[i].x;
    y = line[i].y;
    if((x > 0 && x < SCREEN_WIDTH) && (y > 0 && y < SCREEN_HEIGHT) ){
      PixelShader(surface, line[i], color);
    }
  }
}

void DrawPolygon( const vector<Vertex>& vertices, screen* screen, int t )
{
  int V = vertices.size();
  vector<Pixel> vertexPixels( V );
  for( int i=0; i<V; ++i ){
    VertexShader( vertices[i], vertexPixels[i],t );
  }

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  DrawPolygonRows( leftPixels, rightPixels, screen, triangles[t].color );
}

void PixelsOnLine(Pixel a, Pixel b, vector<Pixel>& line){
  vec2 pix_a = vec2(a.x,a.y);
  vec2 pix_b = vec2(b.x,b.y);
  vec2 delta = glm::abs(pix_a-pix_b);
  int pixels = glm::max(delta.x,delta.y)+1;
  line.resize(pixels);
  Interpolate(a,b,line);
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result){
  int N = result.size();
  vec3 pix_a = vec3(a.x,a.y,a.zinv);
  vec3 pix_b = vec3(b.x,b.y,b.zinv);

  vec4 pixelPosA = a.pos3d;
  vec4 pixelPosB = b.pos3d;

  vec3 step = vec3(pix_b - pix_a) / float(max(N-1, 1));
  vec4 pixel_step = vec4(pixelPosB - pixelPosA) / float(max(N-1,1));
  vec3 current = pix_a;
  vec4 currentPixel = pixelPosA;
  for (int i = 0; i < N; i++){
    result[i].x = round(current.x);
    result[i].y = round(current.y);
    result[i].zinv = current.z;
    result[i].pos3d = currentPixel;

    current += step;
    currentPixel += pixel_step;
  }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels,
vector<Pixel>& leftPixels,vector<Pixel>& rightPixels )
{

// 1. Find max and min y-value of the polygon
// and compute the number of rows it occupies.
  int max_y = -numeric_limits<int>::max();
  int min_y = +numeric_limits<int>::max();

  for (uint32_t i = 0; i < vertexPixels.size(); i++){
    if (vertexPixels[i].y < min_y) {min_y = vertexPixels[i].y;}
    if (vertexPixels[i].y > max_y) {max_y = vertexPixels[i].y;}
  }

  int rows = max_y - min_y + 1;

// 2. Resize leftPixels and rightPixels
// so that they have an element for each row.
  leftPixels.resize(rows);
  rightPixels.resize(rows);

// 3. Initialize the x-coordinates in leftPixels
// to some really large value and the x-coordinates
// in rightPixels to some really small value.

  for(int i = 0; i < rows; i++)  {
    leftPixels[i].x  = +numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();
    leftPixels[i].y = min_y + i;
    rightPixels[i].y = min_y + i;
  }

// 4. Loop through all edges of the polygon and use
// linear interpolation to find the x-coordinate for
// each row it occupies. Update the corresponding
// values in rightPixels and leftPixels.
  // printf("%d \n", vertexPixels.size());
  for (uint32_t i = 0; i < vertexPixels.size(); i++){
    int j = (i+1)%vertexPixels.size();

    vector<Pixel> line(1);
    PixelsOnLine(vertexPixels[i], vertexPixels[j], line);

    //ongoing_y = min_y;
    for (uint32_t j = 0; j < line.size(); j++){
      int offset = line[j].y - min_y;
      if ( leftPixels[offset].x > line[j].x) {
        leftPixels[offset].x = line[j].x;
        leftPixels[offset].zinv = line[j].zinv;
        leftPixels[offset].pos3d = line[j].pos3d;
      }
      if (rightPixels[offset].x < line[j].x) {
        rightPixels[offset].x = line[j].x;
        rightPixels[offset].zinv = line[j].zinv;
        rightPixels[offset].pos3d = line[j].pos3d;
      }

    }
  }


}


int main( int argc, char* argv[] )
{
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  SDL_Renderframe(screen);

  while ( Update())
    {

      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}



vector<Lightning> allocateBolts(vector<Lightning> boltTree){
  int bolt_index;
  int recovery_level = int(0.5f*MAX_BOLTS_ONSCREEN);
  //printf("boltTree size: %d\n", boltTree.size());
  if((boltTree.size() > MAX_BOLTS_ONSCREEN) && boltSwitch == true ){
    boltSwitch = false;
    bolt_index = get_random_index(boltTree.size());
    //printf("Bolt index: %d\n", bolt_index);
  //  boltTree.erase(boltTree.begin()+bolt_index,boltTree.begin()+boltTree.size());
    boltTree.erase(boltTree.begin()+bolt_index);


  }
  else if((boltTree.size() > recovery_level) && boltSwitch == false){
    bolt_index = get_random_index(boltTree.size());
    //printf("Bolt index: %d\n", bolt_index);
    //boltTree.erase(boltTree.begin()+bolt_index,boltTree.begin()+boltTree.size());
    boltTree.erase(boltTree.begin()+bolt_index);

  }
  else{
    boltSwitch = true;
  }
  return boltTree;
}







/*Place your drawing here*/
void Draw(screen* screen)
{
  for( int y=0; y<SCREEN_HEIGHT; ++y ){
    for( int x=0; x<SCREEN_WIDTH; ++x ){
      depthBuffer[y][x] = 0;
    }
  }
/* Clear buffer */
LoadTestModel(triangles);
memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  for( uint32_t i=0; i<triangles.size(); ++i )
  {
    //algorithm that acquires co-ordinates
    //ProduceLightning(boltSegments, vec3(currentNormal));

    vector<Vertex> vertices(3);



    mat4 transformation_mat = RotateAllAxes(pitch,yaw,0.0f);

    currentNormal = triangles[i].normal;
    currentColour = triangles[i].color;

    vertices[0].pos = transformation_mat * vec4(triangles[i].v0.x,triangles[i].v0.y,triangles[i].v0.z,1);
    vertices[0].originalPos = triangles[i].v0;
    vertices[0].triangleIndex = i;
    vertices[0].vertexIndex = 0;

    vertices[1].pos = transformation_mat *vec4(triangles[i].v1.x,triangles[i].v1.y,triangles[i].v1.z,1);
    vertices[1].originalPos = triangles[i].v1;
    vertices[1].triangleIndex = i;
    vertices[1].vertexIndex = 1;

    vertices[2].pos = transformation_mat * vec4(triangles[i].v2.x,triangles[i].v2.y,triangles[i].v2.z,1);
    vertices[2].originalPos = triangles[i].v2;
    vertices[2].triangleIndex = i;
    vertices[2].vertexIndex = 2;

    int triangleIndex = i;
    DrawPolygon( vertices, screen, triangleIndex);
    //if(triangleIndex == 0) printf("Vertex pos: %f %f %f\n", vertices[0].pos.x,vertices[0].pos.y,vertices[0].pos.z);
    // printf("%f %f %f
    //         %f %f %f
    //         %f %f %f", transformation_mat[0][0], transformation_mat[0][1])

  if(boltSwitch){
    vector<Vertex> shockPoints(2);
    vector<Vertex> mirrorPoints(2);
    vector<Vertex> mixBolts(2);
    boltSegments = ProduceLightning(boltSegments,vec3(0.62,-0.1,0.3));
    mirrorBolts =  ProduceLightning(boltSegments,vec3(currentNormal));
    boltAvg = ProduceLightning(boltSegments, 0.5f*(vec3(0.62f+currentNormal.x,
      -0.1f+currentNormal.y, 0.3f+currentNormal.z)));


  //  originBolt = drawOrigin(vec3(0.01,-0.62,-0.6));
    //
    for(uint32_t b = 0; b < boltSegments.size(); b++){
      shockPoints[0].pos = transformation_mat * vec4(boltSegments[b].start.x, boltSegments[b].start.y, boltSegments[b].start.z, 1);
      shockPoints[0].originalPos = vec4(boltSegments[b].start.x, boltSegments[b].start.y, boltSegments[b].start.z, 1);
      shockPoints[1].pos = transformation_mat * vec4(boltSegments[b].end.x, boltSegments[b].end.y, boltSegments[b].end.z, 1);
      shockPoints[0].originalPos = vec4(boltSegments[b].end.x, boltSegments[b].end.y, boltSegments[b].end.z, 1);
      mirrorPoints[0].pos = transformation_mat * vec4(mirrorBolts[b].start.x, mirrorBolts[b].start.y, mirrorBolts[b].start.z, 1);
      mirrorPoints[0].originalPos = vec4(mirrorBolts[b].start.x, mirrorBolts[b].start.y, mirrorBolts[b].start.z, 1);
      mirrorPoints[1].pos = transformation_mat * vec4(mirrorBolts[b].end.x, mirrorBolts[b].end.y, mirrorBolts[b].end.z, 1);
      mirrorPoints[0].originalPos = vec4(mirrorBolts[b].end.x, mirrorBolts[b].end.y, mirrorBolts[b].end.z, 1);
      mixBolts[0].pos = transformation_mat * vec4(boltAvg[b].start.x, boltAvg[b].start.y, boltAvg[b].start.z, 1);
      mixBolts[0].originalPos = vec4(boltAvg[b].start.x, boltAvg[b].start.y, boltAvg[b].start.z, 1);
      mixBolts[1].pos = transformation_mat * vec4(boltAvg[b].end.x, boltAvg[b].end.y, boltAvg[b].end.z, 1);
      mixBolts[0].originalPos = vec4(boltAvg[b].end.x, boltAvg[b].end.y, boltAvg[b].end.z, 1);

      int boltIndex = b;
      currentColour = vec3(1.0f,0.97f,0.01f) * 14.f* indirectIntensity - 5.95f* vec3(0.75f,0.35f,0.15f);
      currentNormal += vec4(0.25f,-0.25f,0.25f,1.0f);
      currentColour *= 0.85f;

      DrawPolygon( shockPoints, screen, boltIndex);
      DrawPolygon( mirrorPoints, screen, boltIndex);
      DrawPolygon( mixBolts, screen, boltIndex);


    }
    // for(uint32_t b = 0; b < originBolt.size(); b++){
    //   soloPoint[0].pos = transformation_mat * vec4(originBolt[b].start.x, originBolt[b].start.y, originBolt[b].start.z, 1);
    //   soloPoint[0].originalPos = vec4(originBolt[b].start.x, originBolt[b].start.y, originBolt[b].start.z, 1);
    //
    //   int boltIndex = b;
    //   currentColour *= 1.5f;
    //   DrawPolygon( soloPoint, screen, boltIndex);
    //
    //
    // }


  }


    ProduceLightning(boltSegments, vec3(currentNormal));

  }
}

vec3 DirectLight(const Pixel &p){
  vec3 pos = vec3(p.pos3d) / p.zinv;
  vec3 diff = vec3(lightPos) - pos;
  float rad = glm::length(diff);
  vec3 normal = vec3(currentNormal);
  vec3 directLight = Intensity(lightPower,rad,diff,normal);
  return directLight;
}

void DrawPolygonRows(const vector<Pixel>& leftPixels,
const vector<Pixel>& rightPixels, screen* s, vec3 col){
    for (uint32_t j = 0; j < leftPixels.size(); j++){
      //Cohen_Suther_Clipping(s, leftPixels[j], rightPixels[j]);
      DrawLineSDL(s, leftPixels[j],rightPixels[j],currentColour);
    }
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
            pitch -= 0.01;
            cameraRot = RotateAllAxes(pitch,0.0f,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2],1);
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2],1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2],1);
            break;
	      case SDLK_DOWN:  /* Move camera backwards */
            pitch += 0.01;
            cameraRot = RotateAllAxes(pitch,0.0f,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2],1);
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2],1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2],1);
            break;
	      case SDLK_LEFT:	 /* Move camera left */
            yaw += 0.01;
            cameraRot = RotateAllAxes(0.0f,yaw,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2],1);
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2],1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2],1);
            break;
        case SDLK_RIGHT: /* Move camera right */
            yaw -= 0.01;
            cameraRot = RotateAllAxes(0.0f,yaw,0.0f);
            right_v   = vec4 (cameraRot[0][0], cameraRot[0][1], cameraRot[0][2],1);
            down_v    = vec4 ( cameraRot[1][0], cameraRot[1][1], cameraRot[1][2],1);
            forward_v = vec4 ( cameraRot[2][0], cameraRot[2][1], cameraRot[2][2],1);
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
          cameraPos = cameraStart;
          lightPos = lightStart;
          pitch = 0;
          yaw = 0;
          right_v = vec4(1,0,0,1);
          down_v = vec4(0,1,0,1);
          forward_v = vec4(0,0,1,1);
          break;
	      }
	  }
    }
  return true;
}



/*
‘std::vector<glm::tvec2<int, (glm::precision)0u> >’ to
‘glm::ivec2 {aka glm::tvec2<int, (glm::precision)0u>}’


*/
