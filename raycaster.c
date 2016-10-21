#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


int line = 1;

int width;
int height;

double cam_width;
double cam_height;


typedef struct {
  int kind; // 0 = plane, 1 = sphere, 2 = camera, 3 = light;
  double diffuse_color[3];
  double specular_color[3];
  //char sym;
  union {
    struct {
      double position[3]; //center position;
      double normal[3];
    } plane;
    struct {
      double position[3]; //center position;
      double radius;
    } sphere;
    struct {
      double color[3];
      double position[3];
      double direction[3];
      double radial_a2;
      double radial_a1;
      double radial_a0;
      double angular_a1;
      double theta;
    } light;
  };
} Object;





typedef struct Pixel {

    unsigned char r, g, b;

 }   Pixel;



typedef struct {
    double width;
    double height;

} Camera;

double clamp(double value){
    if(value > 1){
        return 1;
    }else if(value < 0){
        return 0;
    }else{
        return value;
    }
}
typedef double* V3;

static inline void v3_add(V3 a, V3 b, V3 c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];


}

static inline void v3_subtract(V3 a, V3 b, V3 c) {

  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

static inline void v3_scale(V3 a, double s, V3 c) {

  c[0] = s * a[0];
  c[1] = s * a[1];
  c[2] = s * a[2];



}

static inline double v3_dot(V3 a, V3 b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void v3_cross(V3 a, V3 b, V3 c) {

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];



}

static inline double v3_distance(V3 a, V3 b){
/*
  for(int i = 0; i < 3; i++){
    printf("a value: %f\n", a[i]);
    printf("b value: %f\n", b[i]);
  }
  */
  double dx = b[0] - a[0];
  double dy = b[1] - a[1];
  double dz = b[2] - a[1];
  double c = sqrt(dx*dx + dy*dy + dz*dz);
  //printf("distance value: %d \n", c);
  return c;
}
// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
    ;
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Object** read_scene(char* filename, Object** objects) {

    int c;
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }

  skip_ws(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects
  int index = 0;
  while (1) {
    //printf("Index: %d \n", index);
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      exit(1);
    }
    if (c == '{') {
      skip_ws(json);

      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
	fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
	exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);

      if (strcmp(value, "camera") == 0) {
            if(index ==0 ){
                index = -1;
            }else{
                index--;
            }
      }else if(strcmp(value, "light") == 0){
          objects[index] = malloc(sizeof(Object));
          objects[index]->light.theta = 0;
          objects[index]->light.radial_a2 = 1;
          objects[index]->light.radial_a1 = 0;
          objects[index]->light.radial_a0 = 0;
          objects[index]->light.color[0] = 1;
          objects[index]->light.color[1] = 1;
          objects[index]->light.color[2] = 1;
          objects[index]->kind = (int) 3;
      } else if (strcmp(value, "sphere") == 0) {
          objects[index] = malloc(sizeof(Object));
          objects[index]->diffuse_color[0] = 0;
          objects[index]->diffuse_color[1] = 0;
          objects[index]->diffuse_color[2] = 0;
          objects[index]->specular_color[0] = 1;
          objects[index]->specular_color[1] = 1;
          objects[index]->specular_color[2] = 1;
          objects[index]->kind = (int) 1;
          //printf("Kind: %d \n", objects[index]->kind);
      } else if (strcmp(value, "plane") == 0) {
          objects[index] = malloc(sizeof(Object));
           objects[index]->diffuse_color[0] = 0;
          objects[index]->diffuse_color[1] = 0;
          objects[index]->diffuse_color[2] = 0;
          objects[index]->specular_color[0] = 1;
          objects[index]->specular_color[1] = 1;
          objects[index]->specular_color[2] = 1;
          objects[index]->kind = (int) 0;
          //printf("Kind: %d \n", objects[index]->kind);
      } else {
	fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
	exit(1);
      }

      skip_ws(json);

    while (1) {
	// , }
	c = next_c(json);
	//printf("C at front of loop: %c\n", c);
	if (c == '}') {
	  // stop parsing this object
	  break;
	} else if (c == ',') {
	  // read another field
	  skip_ws(json);
	  char* key = next_string(json);
	  //printf(key);
	  skip_ws(json);
	  expect_c(json, ':');
	  skip_ws(json);
	  if ((strcmp(key, "width") == 0) ||
	      (strcmp(key, "height") == 0) ||
          (strcmp(key, "radial-a0") == 0) ||
          (strcmp(key, "radial-a1") == 0) ||
          (strcmp(key, "radial-a2") == 0) ||
          (strcmp(key, "angular-a1") == 0) ||
          (strcmp(key, "theta") == 0) ||
	      (strcmp(key, "radius") == 0)) {
	    double value = next_number(json);
	    //printf("VALUE: %f", value);
	    //printf("Value: %f\n", value);
	    if(strcmp(key, "width") == 0){
        cam_width = value;
	    }
	    if(strcmp(key, "height") == 0){
        cam_height = value;
	    }
	    if(strcmp(key, "radius") == 0){
          if(objects[index]->kind == 1){
              objects[index]->sphere.radius = value;
              //printf("radius: %f", value);

          }
	    }else if(strcmp(key, "radial-a0") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a0 = value;
                //printf("value: %f \n", value);
            }
	    }else if(strcmp(key, "radial-a1") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a1 = value;
            }
	    }else if(strcmp(key, "radial-a2") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a2 = value;
            }
	    }else if(strcmp(key, "angular-a1") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.angular_a1 = value;
            }
	    }else if(strcmp(key, "theta") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.theta = value;
            }
	    }
	  }else if ((strcmp(key, "diffuse_color") == 0) ||
             (strcmp(key, "specular_color") == 0) ||
             (strcmp(key, "color") == 0) ||
		     (strcmp(key, "position") == 0) ||
             (strcmp(key, "direction") == 0) ||
		     (strcmp(key, "normal") == 0)) {
	    double* value = next_vector(json);
	    if(strcmp(key, "diffuse_color") == 0){
            objects[index]->diffuse_color[0] = value[0];
            objects[index]->diffuse_color[1] = value[1];
            objects[index]->diffuse_color[2] = value[2];
            //printf("got to color...");
         }else if(strcmp(key, "specular_color") == 0){
            objects[index]->specular_color[0] = value[0];
            objects[index]->specular_color[1] = value[1];
            objects[index]->specular_color[2] = value[2];
         }else if(strcmp(key, "color") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.color[0] = value[0];
                objects[index]->light.color[1] = value[1];
                objects[index]->light.color[2] = value[2];
            }
         }else if(strcmp(key, "position") == 0){
            if(objects[index]->kind == 0){
                objects[index]->plane.position[0] = value[0];
                objects[index]->plane.position[1] = value[1]; //*-1
                objects[index]->plane.position[2] = value[2];
                //printf("got to position(plane)...");
            }
            else if(objects[index]->kind == 1){
                objects[index]->sphere.position[0] = value[0];
                objects[index]->sphere.position[1] = value[1]; //*-1
                objects[index]->sphere.position[2] = value[2];
            }else if(objects[index]->kind == 3){
                objects[index]->light.position[0] = value[0];
                objects[index]->light.position[1] = value[1];  //*-1
                objects[index]->light.position[2] = value[2];
            }else{
                //this attribute does not belong in this objects
                printf("THIS WENT WRONG");
            }

         }else if(strcmp(key, "normal") == 0){
            if(objects[index]->kind == 0){
                objects[index]->plane.normal[0] = value[0];
                objects[index]->plane.normal[1] = value[1];
                objects[index]->plane.normal[2] = value[2];
            }else
                //this attribute does not belong in this objects
                printf("This went TERRIBLY wrong...");
            }
         else if (strcmp(key, "direction") == 0){

           if(objects[index]->kind == 3){
                objects[index]->light.direction[0] = value[0];
                objects[index]->light.direction[0] = value[0];
                objects[index]->light.direction[0] = value[0];
           }else{
                //this attribute does not belong in this objects
           }
        }
	   }else {
	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
		    key, line);
	    //char* value = next_string(json);
	  }

	  skip_ws(json);
	  }else {
	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
	  exit(1);
	}
      }
      //printf("C is: %c \n", c);
      skip_ws(json);
      //printf("C is: %c \n", c);
      c = next_c(json);
      //printf("C is: %c \n", c);
      if (c == ',') {
        // noop
        skip_ws(json);
      } else if (c == ']') {
        //printf("Closing File....");
        fclose(json);
        objects[index+1] = NULL;
        //printf("Kind: %d \n", objects[index-2]->kind);
        return objects;
      } else {
	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
	exit(1);
      }
    }
    index++;
  }

}
/////////////////////////////////////////////////////////////////////////////////////


static inline double sqr(double v) {
  return v*v;
}

static inline double to_radians(double degrees){
    return degrees * (M_PI/180.0);
}

static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double sphere_intersection(double* Ro, double* Rd,
			     double* C, double r) {

    double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
    double b = (2*(Rd[0]*(Ro[0] - C[0]) + Rd[1]*(Ro[1] - C[1]) + Rd[2]*(Ro[2] - C[2])));
    double c = sqr(Ro[0] - C[0]) + sqr(Ro[1] - C[1]) + sqr(Ro[2] - C[2]) - sqr(r);

    double det = sqr(b) - 4 *a*c;

    //printf("%f ", det);

  if (det < 0){return -1;}

    //printf("test");
    det = sqrt(det);

  double t0 = (-b - det / (2*a));
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;

}

double plane_intersection(double* Ro, double* Rd, double* L, double* N){
    double D;
    double t;

    for(int i = 0; i <3; i++){
        D += -1*L[i]*N[i];
    }

    t = -(N[0]*Ro[0] + N[1]*Ro[1] + N[2]*Ro[2] + D)/(N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2]);

    return t;

}

Pixel** raycast(Object** scene){
    Pixel** buffer;


    double cx = 0;
    double cy = 0;
    double h = cam_height;
    double w = cam_width;

    int M = width; //width
    int N = height; //height

    width = M;
    height = N;

    buffer = malloc(sizeof(Pixel*)* M );

    double pixheight = h / M;
    double pixwidth = w / N;

    int pixcount = 0;
    for (int y = 0; y < M; y += 1) {
    buffer[y] = malloc(sizeof(Pixel)*N);
    for (int x = 0; x < N; x += 1) {
      double Ro[3] = {0, 0, 0};
      double Rd[3] = {
        (cx - (w/2) + pixwidth * (x + 0.5)),
        (cy - (h/2) + pixheight * (y + 0.5)),
        1
      };
      normalize(Rd);

      double best_t = INFINITY;
      Object* current_object;
      Object* closest_object;
      for (int i=0; scene[i] != 0; i += 1) {
        double t = 0;
        if(scene[i]->kind == 3){
          continue;
        }
    current_object = scene[i];
	switch(current_object->kind) {
	case 1:
	  t = sphere_intersection(Ro, Rd,
				    scene[i]->sphere.position,
				    scene[i]->sphere.radius);
				    //printf("%f", t);
	  break;
    case 0:
        t = plane_intersection(Ro, Rd,
                     scene[i]->plane.position,
                     scene[i]->plane.normal);
        break;
	default:
	  // Horrible error
	  printf("horrible error");
	  exit(1);
	}


	if (t > 0 && t < best_t){
        best_t = t;
        closest_object = current_object;
	}
      }


      //we have the closest object
      //printf("GOT CLOSEST OBJECT");
      double* color = malloc(sizeof(double)*3);
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;

      double* Ron = malloc(sizeof(double)*3);
      double* Rdn = malloc(sizeof(double)*3);

      //shadow test
      for(int i = 0; scene[i] != NULL; i++){ //loop through lights
        //if not a light continue...
        if(scene[i]->kind != 3){
            continue;
        }
        //printf("current light a0: %d \n", scene[i]->light.radial_a0);
        //printf("LOOPING THROUGH LIGHTS \n");

        v3_scale(Rd, best_t, Ron);
        v3_add(Ron, Ro, Ron);
        //printf("GOT RON");
        v3_subtract(scene[i]->light.position, Ron, Rdn);
        //printf("GOT RDN");
        normalize(Rdn);
        double distance_to_light = v3_distance(Ron, scene[i]->light.position);
        //printf("Distance to light: %f\n", distance_to_light);
        //printf("GOT DISTANCE TO LIGHT");

        double new_best_t = INFINITY;
        Object* closest_shadow_object = NULL;
        for(int j = 0; scene[j] != NULL; j++){ //loop through objects
            if(scene[j]->kind == 3){
                continue;
            }
            //printf("LOOPING THROUGH OBJECTS \n");

            if(scene[j] == closest_object){
                //printf("They are equal...");
                continue;
            }
              double new_t = 0;
              Object* current_shadow_object = scene[j];
              switch(scene[j]->kind) {
              case 1:
                new_t = sphere_intersection(Ron, Rdn, scene[j]->sphere.position, scene[j]->sphere.radius);
                break;
              case 0:
                new_t = plane_intersection(Ron, Rdn, scene[j]->plane.position, scene[j]->plane.normal);
                break;
              default:
                // ERROR
                break;
              }


              if (new_t > distance_to_light) {
                continue;
              }

              if (new_t > 0 && new_t < best_t){
                new_best_t = new_t;
                closest_shadow_object = current_shadow_object;
              }




        }

         //
      double frad(double dl, Object* light){
          double a0 = light->light.radial_a0;
          double a1 = light->light.radial_a1;
          double a2 = light->light.radial_a2;

          //printf("A0 = %f \n", a0);
          //printf("A1 = %f \n", a1);
          //printf("A2 = %f \n", a2);
          double ans =  1/(a0 + a1*dl + a2*sqr(dl));
          //printf("returned %f from frad \n", ans);
          return ans;
      }

      double fang(Object* light, double* Vo){ //Vo = N
          //if not spotlight return 1.0
          double a1 = light->light.angular_a1;
          normalize(light->light.direction);
          if(light->light.theta == 0){
            //printf("returned 1 in fang \n");
            return 1;
          }else if(v3_dot(Vo, light->light.direction) < cos(to_radians(light->light.theta))){ //not in spotlight
            //printf("returned 0 in fang \n");
            return 0;
          }else{
            //printf("returned other value in fang \n");
            return pow(v3_dot(Vo, light->light.direction), a1);
          }

      }

      double diffuse(double light_color, double object_color, double* N, double* L){

            double n_dot_l = v3_dot(N, L);
            double result;
            //printf("N dot L is: %f\n", n_dot_l);
            if(n_dot_l > 0){

                //printf("light color = %f \n", light_color);
                //printf("object folor = %f \n", object_color);
                result =  object_color*light_color*n_dot_l;
                //printf("DIFFUSE VALUE: %f \n", result);
                return result;
            }else{
                return 0;
            }

      }

      double specular(double light_color, double specular_color, double* V, double* R, double* N, double* L){

            double v_dot_r = v3_dot(V, R);
            double n_dot_l = v3_dot(N, L);
            double ns = 100;
            double result;

            if(v_dot_r > 0 && n_dot_l > 0){
                //printf("V DOT R: %f \n", v_dot_r);
                //printf("N DOT L: %f \n", n_dot_l);
                result = specular_color*pow(v_dot_r, ns);
                if(result > 0){
                    //printf("Specular value: %f \n", result);
                }

                return result;
            }else{
                //printf("zero spec value...\n");
                return 0;
                }
      }
      if (closest_shadow_object == NULL) { //not in shadow
        //printf("distance to light: %d\n", distance_to_light);
        //printf("NOT IN SHADOW \n");
        // N, L, R, V
        double* Vobject = malloc(sizeof(double)*3);
        v3_subtract(Ron, scene[i]->light.position, Vobject);
        v3_scale(Vobject, -1, Vobject);
        normalize(Vobject);
        double* N = malloc(sizeof(double)*3);
        double* L = malloc(sizeof(double)*3);
        double* R = malloc(sizeof(double)*3);
        double* V = malloc(sizeof(double)*3);
        if(closest_object->kind == 0){
            //printf("THIS IS A PLANE \n");
            N[0] = closest_object->plane.normal[0];
            N[1] = closest_object->plane.normal[1];
            N[2] = closest_object->plane.normal[2];
            v3_scale(N, -1, N);  // plane
            normalize(N);
        }else if(closest_object->kind == 1){
            //printf("THIS IS A SPHERE \n");
            v3_subtract(Ron, closest_object->sphere.position, N); // sphere
            normalize(N);
        }else{
            printf("error...");
            //something went wrong somewhere :/
        }
        v3_scale(Rdn, -1, L); // light_position - Ron;
        //normalize(L);

        v3_scale(N, 2, R);
        v3_scale(N, v3_dot(R, L), R);
        v3_subtract(R, L, R);//reflection of L

        V[0] = Rd[0];
        V[1] = Rd[1];
        V[2] = Rd[2];
        //v3_scale(Rd, -1, V);
        //normalize(V);
        //normalize(R);


         // uses object's diffuse color
         // uses object's specular color


        color[0] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[0], closest_object->diffuse_color[0], N, L) + specular(scene[i]->light.color[0], closest_object->specular_color[0], V, R, N, L));
        color[1] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[1], closest_object->diffuse_color[1], N, L) + specular(scene[i]->light.color[1], closest_object->specular_color[1], V, R, N, L));
        color[2] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[2], closest_object->diffuse_color[2], N, L) + specular(scene[i]->light.color[2], closest_object->specular_color[2], V, R, N, L));



        free(N);
        free(L);
        free(R);
        free(V);

      }

      }



      int temp = y;
      buffer[temp][x].r = clamp(color[0])*255;
      buffer[temp][x].g = clamp(color[1])*255;
      buffer[temp][x].b = clamp(color[2])*255;




/*
      if (best_t > 0 && best_t != INFINITY) {
            //printf("%c", best_color); //hit
            buffer[pixcount].r = best_color[0]*255;
            buffer[pixcount].g = best_color[1]*255;
            buffer[pixcount].b = best_color[2]*255;
      } else {
            buffer[pixcount].r = 0;
            buffer[pixcount].g = 0;
            buffer[pixcount].b = 0;
      }
*/
        pixcount++;
        //printf("%d", pixcount);
    }
    //printf("\n");
  }
  for(int y = 0; y < height; y++){
    for(int x = 0; x < width; x++){
        Pixel temp;
        temp.r = buffer[height-y-1][x].r;
        temp.g = buffer[height-y-1][x].g;
        temp.b = buffer[height-y-1][x].b;
        buffer[height-y-1][x] = buffer[y][x];
        buffer[y][x].r = temp.r;
        buffer[y][x].g = temp.g;
        buffer[y][x].b = temp.b;
    }
  }
  return buffer;
}


int writeP6(char* fname, Pixel** buffer){
    FILE* fh;
    fh = fopen(fname, "wb");

    fprintf(fh, "P6\n");
    fprintf(fh, "%d ", width);
    fprintf(fh, "%d\n", height);
    fprintf(fh, "%d\n", 255);


    for(int i = 0; i < width; i++){
        for(int j = 0; j<height; j++){
        unsigned char* rgb;
        rgb = malloc(sizeof(unsigned char)*64);
        rgb[0] = buffer[i][j].r;
        rgb[1] = buffer[i][j].g;
        rgb[2] = buffer[i][j].b;
        fwrite(rgb, 1, 3, fh);
        }


    };

    fclose(fh);


}

int main(int argc, char *argv[]){
    //printf("%d", argc);

    if(argc != 5){
        fprintf(stderr, "Incorrect usage of arguments: width height input.json output.ppm");
        return 1;
    }


    width = atoi(argv[1]);

    height = atoi(argv[2]);

    Object** objects;
    objects = malloc(sizeof(Object*)*144);

    Object** scene = read_scene(argv[3], objects);

    Pixel** buffer = raycast(scene);

    writeP6(argv[4], buffer);


  return 0;

}
