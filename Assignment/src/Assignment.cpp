#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <utility>
#include <fstream>
#include <vector>
#include <functional>
#include <glm/glm.hpp>
// #include <glm/gtx/hash.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>
#include <ctime>

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480
#define LIGHT 45
#define REFRACTIVE_INDEX 1.3

#define pi 3.14159265359

RayTriangleIntersection get_closest_reflection(vec3 int_point, vec3 direction, vector<ModelTriangle> &triangles, int index, int depth);
RayTriangleIntersection get_closest_refraction(vec3 int_point, vec3 direction, vector<ModelTriangle> &triangles, int index, int depth);

vec3     o(0.0, 0.0, 0.0); // origin
// vec3   cam(140.0,135.0,210.0); // cornell cam
vec3   cam(0.0, 0.0, 4.0); // sphere cam
vector<vec3> lights{
	// vec3(-0.1, 0.8, -0.1),
	// vec3(-0.1, 0.8, 0.0),
	// vec3(-0.1, 0.8, 0.1),
	// vec3(0.0, 0.8, -0.1),
	// vec3(0.0, 0.8, 0.0),
	// vec3(0.0, 0.8, 0.1),
	// vec3(0.1, 0.8, -0.1),
	// vec3(0.1, 0.8, 0.0),
	// vec3(0.1, 0.8, 0.1),
	// vec3(-0.1, 0.9, -0.1),
	// vec3(-0.1, 0.9, 0.0),
	// vec3(-0.1, 0.9, 0.1),
	// vec3(0.0, 0.9, -0.1),
	// vec3(0.0, 0.9, 0.0),
	// vec3(0.0, 0.9, 0.1),
	// vec3(0.1, 0.9, -0.1),
	// vec3(0.1, 0.9, 0.0),
	// vec3(0.1, 0.9, 0.1),
	vec3(-0.2, 1.0, 1.8),
	vec3(-0.2, 1.0, 2.0),
	vec3(-0.2, 1.0, 2.2),
	vec3(0.0, 1.0, 1.8),
	vec3(0.0, 1.0, 2.0),
	vec3(0.0, 1.0, 2.2),
	vec3(0.2, 1.0, 1.8),
	vec3(0.2, 1.0, 2.0),
	vec3(0.2, 1.0, 2.2)
};
// std::vector<glm::vec3> lights = {
// 	glm::vec3(0.0,0.85,0.0),
// 	glm::vec3(0.0,0.85,0.1),
// 	glm::vec3(0.0,0.85,-0.1),
// 	glm::vec3(0.1,0.85,0.0),
// 	glm::vec3(-0.1,0.85,0.1),
// 	glm::vec3(0.1,0.85,-0.1),
// 	glm::vec3(0.1,0.85,0.1),
// 	glm::vec3(-0.1,0.85,0.0),
// 	glm::vec3(-0.1,0.85,-0.1)
// };
vec3 *light_draw = &lights[4]; // 4 23 cornell light - 1.0 is proportional to scale used when loading obj
// vec3 light(0.0, 200, 00); // sphere light - values are proportional to scale used when loading obj
mat3 cam_orientation(vec3(1.0,0.0,0.0),vec3(0.0,1.0,0.0),vec3(0.0,0.0,1.0));

float focal = 500.0;
bool orbiting = false, show_light = false;
// bool proximity = true, angle_of = true, shadows = true, specular = true; // cornell values
bool proximity = true, angle_of = true, shadows = true, specular = true; // sphere values

void draw(DrawingWindow &window) {
	window.clearPixels();
	// for (size_t y = 0; y < window.height; y++) {
	// 	for (size_t x = 0; x < window.width; x++) {
	// 		float red = rand() % 256;
	// 		float green = 0.0;
	// 		float blue = 0.0;
	// 		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	// 		// window.setPixelColour(x, y, colour);
	// 	}
	// }
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

vector<CanvasPoint> interpolate_points(CanvasPoint from, CanvasPoint to, int steps){
	vector<CanvasPoint> points;
	float x_step = (to.x - from.x)/(steps-1);
	float y_step = (to.y - from.y)/(steps-1);
	float d_step = (to.depth - from.depth)/(steps-1);

	CanvasPoint temp = from;
	points.push_back(temp);

	for (int i = 0; i < steps-1; i++) {
		temp.x = temp.x + x_step;
		temp.y = temp.y + y_step;
		temp.depth = temp.depth + d_step;
		points.push_back(temp);
	}
	return points;
}

vector<TexturePoint> interpolate_points(TexturePoint from, TexturePoint to, int steps){
	vector<TexturePoint> points;
	float x_step = (to.x - from.x)/(steps-1);
	float y_step = (to.y - from.y)/(steps-1);

	TexturePoint temp = from;
	points.push_back(temp);

	for (int i = 0; i < steps-1; i++) {
		temp.x = temp.x + x_step;
		temp.y = temp.y + y_step;
		points.push_back(temp);
	}
	return points;
}

void draw_line(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window) {
	float x_diff = to.x - from.x;
	float y_diff = to.y - from.y;
	float steps = std::max(abs(x_diff),abs(y_diff));

	float x_step = x_diff/steps;
	float y_step = y_diff/steps;

	uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

	for(float i = 0.0; i < steps; i++) {
		int x = round(from.x + (x_step * i));
		int y = round(from.y + (y_step * i));
		if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
			window.setPixelColour(x, y, c);
		}
	}
}

void draw_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
	CanvasPoint a = triangle[0];
	CanvasPoint b = triangle[1];
	CanvasPoint c = triangle[2];

	draw_line(a, b, colour, window);
	draw_line(b, c, colour, window);
	draw_line(c, a, colour, window);
}

CanvasPoint find_mid(CanvasPoint top, CanvasPoint mid, CanvasPoint bot) {
	float x_mid = top.x+((mid.y-top.y)/(bot.y-top.y))*(bot.x-top.x);
	return CanvasPoint(round(x_mid), mid.y);
}

void fill_half_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, vector<vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];

	vector<CanvasPoint> left = interpolate_points(top, mid, abs(mid.y-top.y)+2);
	vector<CanvasPoint> right = interpolate_points(top, bot, abs(mid.y-top.y)+2);

	for(int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);

		vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+2);
		for(int j = 0; j < points.size(); j++) {
			int x = round(points[j].x);
			int y = round(points[j].y);
			if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
				if(-1/points[j].depth > depths[x][y]) {
					depths[x][y] = -1/points[j].depth;
					uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
					window.setPixelColour(x,y,c);
				}
			}
		}
	}
}

void texture_half_triangle(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window, vector<vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];
	
	vector<CanvasPoint> left = interpolate_points(top, mid, abs(mid.y-top.y)+2);
	vector<CanvasPoint> right = interpolate_points(top, bot, abs(mid.y-top.y)+2);
	vector<TexturePoint> left_texture = interpolate_points(top.texturePoint, mid.texturePoint, abs(mid.y-top.y)+2);
	vector<TexturePoint> right_texture = interpolate_points(top.texturePoint, bot.texturePoint, abs(mid.y-top.y)+2);

	for(int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);

		vector<CanvasPoint> points = interpolate_points(left[i], right[i], steps+2);
		vector<TexturePoint> points_texture = interpolate_points(left_texture[i], right_texture[i], steps+2);
		for(int j = 0; j < points.size(); j++) {
			int x = round(points[j].x);
			int y = round(points[j].y);
			if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
				if(-1/points[j].depth > depths[x][y]) {
					depths[x][y] = -1/points[j].depth;
					window.setPixelColour(x, y, texture.pixels[round(points_texture[j].y)*texture.width + round(points_texture[j].x)]);
				}
			}
		}
	}
}

void fill_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window, vector<vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];

    if (bot.y < mid.y) {
        std::swap(bot, mid);
    } if (mid.y < top.y) {
        std::swap(mid, top);
    } if (bot.y < mid.y) {
        std::swap(bot, mid);
    }

	CanvasPoint mid_2 = find_mid(top, mid, bot);
	mid_2.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	CanvasTriangle t_1 = CanvasTriangle(top,mid,mid_2);
	CanvasTriangle t_2 = CanvasTriangle(bot,mid,mid_2);

	fill_half_triangle(t_1, colour, window, depths);
	fill_half_triangle(t_2, colour, window, depths);

	// draw_triangle(triangle, Colour(255,255,255), window);
}

void texture_triangle(TextureMap texture, CanvasTriangle triangle, DrawingWindow &window, vector<vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];

    if (bot.y < mid.y) {
        std::swap(bot, mid);
    } if (mid.y < top.y) {
        std::swap(mid, top);
    } if (bot.y < mid.y) {
        std::swap(bot, mid);
    }

	CanvasPoint mid_2 = find_mid(top, mid, bot);
	mid_2.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	float scale = (mid.y - top.y)/(bot.y-top.y);
	mid_2.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	mid_2.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	CanvasTriangle t_1 = CanvasTriangle(top,mid,mid_2);
	CanvasTriangle t_2 = CanvasTriangle(bot,mid,mid_2);

	texture_half_triangle(t_1, texture, window, depths);
	texture_half_triangle(t_2, texture, window, depths);
	
	// draw_triangle(triangle, Colour(255,255,255), window);
}

void draw_wireframe(vector<ModelTriangle> &triangles, DrawingWindow &window, TextureMap &texture) {

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle triangle = triangles[i];
		CanvasTriangle t;
		for(int j = 0; j < triangles[i].vertices.size(); j++) {
			vec3 vertex = triangles[i].vertices[j];
			vec3 cam_to_vertex(vertex.x - cam.x, vertex.y - cam.y, vertex.z - cam.z);
			vec3 adjusted_vertex = cam_to_vertex * cam_orientation;
		
			int u = -(focal * (adjusted_vertex.x)/(adjusted_vertex.z)) + (window.width/2);
			int v = (focal * (adjusted_vertex.y)/(adjusted_vertex.z)) + (window.height/2);
			;
			t.vertices[j] = CanvasPoint(u,v, adjusted_vertex.z);
			t.vertices[j].texturePoint = triangles[i].texturePoints[j];
		}
		draw_triangle(t, Colour(255,255,255), window);
	}

	vec3 cam_to_vertex = vec3(light_draw->x - cam.x, light_draw->y - cam.y, light_draw->z - cam.z);
    vec3 adjusted_vector = cam_to_vertex * cam_orientation;

    int u = -(focal * (adjusted_vector.x)/(adjusted_vector.z)) + (window.width / 2);
    int v = (focal * (adjusted_vector.y)/(adjusted_vector.z)) + (window.height / 2);

    // prints red pixels to show light location
	if(u > 0 && u < window.width-1 && v > 0 && v < window.height-1) {
		window.setPixelColour(u,   v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
		window.setPixelColour(u+1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
		window.setPixelColour(u, v+1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
		window.setPixelColour(u-1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
		window.setPixelColour(u, v-1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	}
}

void draw_rasterise(vector<ModelTriangle> &triangles, DrawingWindow &window, TextureMap &texture) {

	vector<vector<float>> depths(window.width, vector<float>(window.height, -(numeric_limits<float>::infinity())));

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle triangle = triangles[i];
		CanvasTriangle t;
		for(int j = 0; j < triangles[i].vertices.size(); j++) {
			vec3 vertex = triangles[i].vertices[j];
			vec3 cam_to_vertex(vertex.x - cam.x, vertex.y - cam.y, vertex.z - cam.z);
			vec3 adjusted_vertex = cam_to_vertex * cam_orientation;
		
			int u = -(focal * (adjusted_vertex.x)/(adjusted_vertex.z)) + (window.width/2);
			int v = (focal * (adjusted_vertex.y)/(adjusted_vertex.z)) + (window.height/2);
			;
			t.vertices[j] = CanvasPoint(u,v, adjusted_vertex.z);
			t.vertices[j].texturePoint = triangles[i].texturePoints[j];
		}
		if (triangles[i].colour.name != "") {
			TextureMap texturex("logo.ppm");//triangle.colour.name);
			for(int j = 0; j < t.vertices.size(); j++) {
				t.vertices[j].texturePoint.x *= texture.width;
				t.vertices[j].texturePoint.y *= texture.height;
			}
			cout << texture.pixels.size() << endl;
			texture_triangle(texturex, t, window, depths);
		} else {
			fill_triangle(t, triangles[i].colour, window, depths);
		}
	}
}

pair<bool,bool> is_shadow(RayTriangleIntersection intersect, vec3 &light, vector<ModelTriangle> &triangles) {

	vec3 shadow_ray = light - intersect.intersectionPoint;

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle tri = triangles[i];

		vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		vec3 sp_vector = intersect.intersectionPoint - triangles[i].vertices[0];
		mat3 de_matrix(-normalize(shadow_ray), e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(t < glm::length(shadow_ray) && t > 0.01f && i != intersect.triangleIndex) {
				bool refract = (triangles[i].refract) ? true : false;
				return std::make_pair(true,refract);//true;
			}
		}
	}
	return std::make_pair(false,false);//false;
}

float get_scale(RayTriangleIntersection rt_int, vec3 &light, int scale) {

	vec3 normal = normalize(rt_int.intersectedTriangle.normal);
	vec3 light_ray = light - rt_int.intersectionPoint;
	vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
	view_ray = normalize(cam_orientation * view_ray);
	vec3 reflection_ray = normalize(normalize(light_ray) - (normal * 2.0f * dot(normalize(light_ray), normal)));

	float scale_p = (proximity) ? LIGHT/(4*pi*(pow(length(light_ray),2))) : 0;
	float scale_a = dot(normal, normalize(light_ray));
	float scale_s = pow(dot(reflection_ray, view_ray),scale);

	if(scale_a > 0 && angle_of) scale_p *= scale_a;
	if(scale_s > 0 && specular) scale_p += scale_s;
	return (scale_p < 1) ? scale_p : 1;
}

float gourad(RayTriangleIntersection rt_int, vec3 &light, int scale) {

	ModelTriangle t = rt_int.intersectedTriangle;
	vec3 light_ray = light - rt_int.intersectionPoint;
	vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
	view_ray = normalize(cam_orientation * view_ray);

	vector<float> scales;
	vector<vec3> reflections;
	for(int i = 0; i < t.normals.size(); i++) {
		vec3 reflection_ray = normalize(normalize(light_ray) - (t.normals[i] * 2.0f * dot(normalize(light_ray), t.normals[i])));
		reflections.push_back(reflection_ray);

		float temp_a = (angle_of) ? dot(t.normals[i], normalize(light_ray)) : 1;
		float temp_p = (proximity) ? LIGHT*temp_a/(4*pi*(pow(length(light_ray),2))) : 0;
		float temp_s = pow(dot(normalize(reflection_ray), view_ray),scale);
		if(temp_s > 0 && specular) temp_p += temp_s;
		scales.push_back(temp_p);
	}
	// vec3 reflection_ray = (1 - rt_int.u - rt_int.v) * reflections[0] + rt_int.u * reflections[1] + rt_int.v * reflections[2];

	float bright = (1-rt_int.u-rt_int.v) * scales[0] + rt_int.u * scales[1] + rt_int.v * scales[2];

	// float scale_p = (proximity) ? LIGHT*scale_a/(4*pi*(pow(length(light_ray),2))) : 0;
	// if(scale_a > 0 && angle_of) scale_p *= scale_a;
	// if(scale_s > 0 && specular) scale_p += (scale_s * 0.2);

	return (bright < 1) ? bright : 1;
}

float phong(RayTriangleIntersection rt_int, vec3 &light, int scale) {

	ModelTriangle t = rt_int.intersectedTriangle;
	vec3 normal = (1 - rt_int.u - rt_int.v) * t.normals[0] + rt_int.u * t.normals[1] + rt_int.v * t.normals[2];
	vec3 light_ray = light - rt_int.intersectionPoint;
	vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
	view_ray = normalize(cam_orientation * view_ray);
	vec3 reflection_ray = normalize(normalize(light_ray) - (normal * 2.0f * dot(normalize(light_ray), normal)));

	float scale_a = (angle_of) ? dot(normal, normalize(light_ray)) : 1;
	float scale_p = (proximity) ? LIGHT*scale_a/(4*pi*(pow(length(light_ray),2))) : 0;
	float scale_s = pow(dot(reflection_ray, view_ray),scale);

	// if(scale_a > 0 && angle_of) scale_p *= scale_a;
	if(scale_s > 0 && specular) scale_p += (scale_s);
	return (scale_p < 1) ? scale_p : 1;
}

uint32_t get_texture(RayTriangleIntersection rt_int, TextureMap &texture) {
	ModelTriangle t = rt_int.intersectedTriangle;

	float x = ((1 - rt_int.u - rt_int.v) * t.texturePoints[0].x + rt_int.u * t.texturePoints[1].x + rt_int.v * t.texturePoints[2].x);
	float y = ((1 - rt_int.u - rt_int.v) * t.texturePoints[0].y + rt_int.u * t.texturePoints[1].y + rt_int.v * t.texturePoints[2].y);

	x *= texture.width;
	y *= texture.height;
	// cout << x << "," << y << endl;
	return texture.pixels[round(y)*texture.width + round(x)];
}

vec3 refract(vec3 incidence, vec3 n, float index) {

	// index = 2.0f - index;
	float cosi = dot(n,incidence);
	// vec3 refracted = (incidence * index - n * (-cosi + index * cosi));
	// return normalize(refracted);

	float eta = index;
	if(cosi > 0) {
		eta = 1.0/eta;
	}
	vec3 refracted = incidence * eta - n * (-cosi + eta * cosi);
	return normalize(refracted);

	// vec3 normal = n;
    // float d = dot(incidence, n);
    // float refr_1 = 1;
    // float refr_2 = index;
    // if (d < 0.0f) {
    //     // outside surface
    //     d = -d;
    // } else {
    //     // inside surface
    //     swap(refr_2, refr_1);
    //     n = -n;
    // }
    // float refr = refr_1 / refr_2;
	// float k = 1 - refr * refr * (1 - d * d);
	// if(k < 0) {
	// 	vec3 reflected_ray = incidence - (2.0f*dot(incidence, normal)*normal);
	// 	return normalize(reflected_ray);
	// 	// cout << "totalinternalreflection" << endl;
	// 	// return vec3(0,0,0);//normalize(incidence - (n * 2.0f * d));
	// }
    // vec3 refracted_ray = normalize(refr * incidence + (refr * d - sqrtf(k)) * n);
    // vec3 refracted_ray = refr * incidence - (refr * d) * n;
    // return refracted_ray;
}
float fresnel(const vec3 I, const vec3 N, const float ior) { 
    float etai = 1, etat = ior; 
    float cosi = dot(I, N);
	if (cosi < 0) cosi = -cosi;
    else swap(etai, etat);
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi)); 
    // Total internal reflection
    if (sint >= 1) { 
        return 1.0; 
    } 
    else { 
        float cost = sqrtf(std::max(0.f, 1 - sint * sint)); 
        cosi = fabsf(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        return (Rs * Rs + Rp * Rp) / 2; 
    } 
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}
RayTriangleIntersection get_closest_reflection(vec3 int_point, vec3 direction, vector<ModelTriangle> &triangles, int index, int depth) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle tri = triangles[i];
		
		vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		vec3 sp_vector = int_point - triangles[i].vertices[0];
		mat3 de_matrix(-direction, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(rti.distanceFromCamera > t && t > 0.0001f && i != index) {
				rti.distanceFromCamera = t;
				rti.intersectedTriangle = triangles[i];
				// rti.intersectedTriangle.normal = cross(e1,e0);
				rti.triangleIndex = i;
				rti.u = u;
				rti.v = v;

				vec3 intersect = triangles[i].vertices[0]+u*e0+v*e1;
				rti.intersectionPoint = intersect;
			}
		}
	}
	if(depth > 4) {
		if(rti.distanceFromCamera )
		rti.intersectedTriangle.colour = Colour(255,255,255);
		return rti;
	}
	else if(rti.intersectedTriangle.mirror) {
		vec3 normal = normalize(rti.intersectedTriangle.normal);
		vec3 reflection_ray = normalize(direction - (normal * 2.0f * dot(direction, normal)));
		rti = get_closest_reflection(rti.intersectionPoint, reflection_ray, triangles, rti.triangleIndex, depth++);
	} 
	else if(rti.intersectedTriangle.refract) {
		vec3 normal = normalize(rti.intersectedTriangle.normal);
		vec3 refracted_ray = refract(direction, normal, REFRACTIVE_INDEX);
		// if(refracted_ray == vec3(0,0,0)) {
		// 	vec3 reflected_ray = direction - (2.0f*dot(direction, normal)*normal);
		// 	reflected_ray = normalize(reflected_ray);
		// 	RayTriangleIntersection reflection = get_closest_refraction(rti.intersectionPoint, refracted_ray, triangles, rti.triangleIndex);
		// 	rti = reflection;
		// } else {
			rti = get_closest_refraction(rti.intersectionPoint, refracted_ray, triangles, rti.triangleIndex, depth++);
		// }
	}
	return rti;
}
RayTriangleIntersection get_closest_refraction(vec3 int_point, vec3 direction, vector<ModelTriangle> &triangles, int index, int depth) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle tri = triangles[i];
		
		vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		vec3 sp_vector = int_point - triangles[i].vertices[0];
		mat3 de_matrix(-direction, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(rti.distanceFromCamera > t && t > 0.0001f && i != index) {
				rti.distanceFromCamera = t;
				rti.intersectedTriangle = triangles[i];
				// rti.intersectedTriangle.normal = cross(e1,e0);
				rti.triangleIndex = i;
				rti.u = u;
				rti.v = v;

				vec3 intersect = triangles[i].vertices[0]+u*e0+v*e1;
				rti.intersectionPoint = intersect;
			}
		}
	}
	if(depth > 4) {
		rti.intersectedTriangle.colour = Colour(255,255,255);
		return rti;
	}
	else if(rti.intersectedTriangle.refract) {
		vec3 normal = normalize(rti.intersectedTriangle.normal);
		vec3 refracted_ray = refract(direction, normal, REFRACTIVE_INDEX);
		// if(refracted_ray == vec3(0,0,0)) {
		// 	vec3 reflected_ray = direction - (2.0f*dot(direction, normal)*normal);
		// 	reflected_ray = normalize(reflected_ray);
		// 	RayTriangleIntersection reflection = get_closest_refraction(rti.intersectionPoint, refracted_ray, triangles, rti.triangleIndex);
		// 	rti = reflection;
		// } else {
			rti = get_closest_refraction(rti.intersectionPoint, refracted_ray, triangles, rti.triangleIndex, depth++);
		// }
	}
	else if(rti.intersectedTriangle.mirror) {
		vec3 normal = normalize(rti.intersectedTriangle.normal);
		vec3 reflection_ray = normalize(direction - (normal * 2.0f * dot(direction, normal)));
		rti = get_closest_reflection(rti.intersectionPoint, reflection_ray, triangles, rti.triangleIndex, depth++);
	} 
	return rti;
}
RayTriangleIntersection get_closest_intersection(vec3 direction, vector<ModelTriangle> &triangles, TextureMap &texture) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = numeric_limits<float>::infinity();
	float dist = numeric_limits<float>::infinity();
	vec3 ray = cam - direction;
	ray = normalize(cam_orientation * ray);

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle tri = triangles[i];
		
		vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
		vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
		vec3 sp_vector = cam - triangles[i].vertices[0];
		mat3 de_matrix(-ray, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(dist > t && t > 0) {
				vec3 intersect = triangles[i].vertices[0]+u*e0+v*e1;
				if(triangles[i].mirror) {
					vec3 normal = normalize(triangles[i].normal);
					vec3 reflection_ray = normalize(ray - (normal * 2.0f * dot(ray, normal)));
					
					RayTriangleIntersection t_rti = get_closest_reflection(intersect, reflection_ray, triangles, i, 1);
					rti = t_rti;
					// if(isinf(t_rti.distanceFromCamera)) rti.inf = true;
				} else if(triangles[i].refract) {
					vec3 normal = normalize(triangles[i].normal);
					vec3 reflection_ray = normalize(ray - (normal * 2.0f * dot(ray, normal)));
					vec3 refracted_ray = refract(ray, normal, REFRACTIVE_INDEX);
					// if(refracted_ray == vec3(0,0,0)) {
						// vec3 reflected_ray = ray - (2.0f*dot(ray, normal)*normal);
						// reflected_ray = normalize(reflected_ray);
						// RayTriangleIntersection reflection = get_closest_refraction(intersect, reflection_ray, triangles, i);
						// rti = reflection;
					// } else {
					RayTriangleIntersection t_rti = get_closest_refraction(intersect, refracted_ray, triangles, i, 1);
					// }
					RayTriangleIntersection r_rti = get_closest_reflection(intersect, reflection_ray, triangles, i, 1);
					Colour refract = t_rti.intersectedTriangle.colour;
					Colour reflect = r_rti.intersectedTriangle.colour;
					if(t_rti.intersectedTriangle.colour.name != "") {
						uint32_t t = get_texture(t_rti, texture);
						float r = (t >> 16) & 0xff;
						float g = (t >>  8) & 0xff;
						float b =         t & 0xff;
						refract = Colour(int(r), int(g), int(b));
					}
					if(r_rti.intersectedTriangle.colour.name != "") {
						uint32_t t = get_texture(r_rti, texture);
						float r = (t >> 16) & 0xff;
						float g = (t >>  8) & 0xff;
						float b =         t & 0xff;
						reflect = Colour(int(r), int(g), int(b));
					}
					float kr = fresnel(ray, normal, REFRACTIVE_INDEX);
					float rfr = 1-kr, rfl = kr;//1-dot(ray, normal), rfl = dot(ray, normal);
					Colour avg(int((refract.red*rfr)+(reflect.red*rfl)),int((refract.green*rfr)+(reflect.green*rfl)),int((refract.blue*rfr)+(reflect.blue*rfl)));
					// t_rti.intersectionPoint = intersect;
					t_rti.intersectedTriangle.colour = avg;
					rti = t_rti;

				} else {
					rti.triangleIndex = i;
					rti.u = u;
					rti.v = v;

					rti.distanceFromCamera = t;
					rti.intersectedTriangle = triangles[i];
					// rti.intersectedTriangle.normal = cross(e1,e0);
					
					rti.intersectionPoint = intersect;
				}
				dist = t;
			}
		}
	}
	return rti;
}

function<float(RayTriangleIntersection rt_int, vec3 &light, int scale)> brightness = phong;

void draw_raytrace(vector<ModelTriangle> &triangles, DrawingWindow &window, TextureMap &texture) {
	// TextureMap texture("logo.ppm");
	for(int x = 0; x < window.width; x++) {
		for(int y = 0; y < window.height; y++) {
			RayTriangleIntersection rt_int = get_closest_intersection(vec3((int(window.width)/2)-x,y-(int(window.height)/2), focal), triangles, texture);

			float scale = 0;
			for(int i = 0; i < lights.size(); i++) {
				std::pair<bool,bool> shadow = is_shadow(rt_int, lights[i], triangles);
				// if(shadow.second) cout << "lighter shadow" << endl;
				float scale_s = (shadow.second) ? 0.1 : 0.1;
				scale += (shadow.first && shadows) ? scale_s : brightness(rt_int, lights[i], 64);
				// scale += is_shadow(rt_int, lights[i], triangles, brightness(rt_int, lights[i], 64));
			}
			scale /= lights.size();
			scale = (scale > 0.1) ? scale : 0.1;
			// if(x%100 == 0 && y%100 == 0) cout << scale << endl << scale_s << endl << endl;
			if(!isinf(rt_int.distanceFromCamera)){
				Colour colour = rt_int.intersectedTriangle.colour;
				bool has_texture = (rt_int.intersectedTriangle.colour.name != "");
				uint32_t c = (255 << 24) + (int(colour.red*scale) << 16) + (int(colour.green*scale) << 8) + int(colour.blue*scale);
				if(has_texture) {
					uint32_t t = get_texture(rt_int, texture);
					float r = (t >> 16) & 0xff;
					float g = (t >>  8) & 0xff;
					float b =         t & 0xff;

					c = (255 << 24) + (int(r*scale) << 16) + (int(g*scale) << 8) + int(b*scale);
				}
				if(rt_int.inf) {
					uint32_t b = (255 << 24) + (0 << 16) + (0 << 8) + 0;
					window.setPixelColour(x,y,b);
				} else window.setPixelColour(x,y,c);
			}
		}
	}
	// if(show_light) {
	// 	vec3 cam_to_vertex = vec3(light_draw->x - cam.x, light_draw->y - cam.y, light_draw->z - cam.z);
	// 	vec3 adjusted_vector = cam_to_vertex * cam_orientation;

	// 	int u = -(focal * (adjusted_vector.x)/(adjusted_vector.z)) + (window.width / 2);
	// 	int v = (focal * (adjusted_vector.y)/(adjusted_vector.z)) + (window.height / 2);

	// 	window.setPixelColour(u,   v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	// 	window.setPixelColour(u+1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	// 	window.setPixelColour(u, v+1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	// 	window.setPixelColour(u-1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	// 	window.setPixelColour(u, v-1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	// }
}

void vertex_normals(vector<ModelTriangle> &triangles) {

	for(int i = 0; i < triangles.size(); i++) {
		// ModelTriangle t = triangles[i];
		vector<vec3> normals;
		for(int v = 0; v < triangles[i].vertices.size(); v++) {
			vec3 vertex = triangles[i].normal;
			int count = 1;
			for(int j = 0; j < triangles.size(); j++) {
				// ModelTriangle t_ = triangles[j];
				for(int u = 0; u < triangles[j].vertices.size(); u++) {
					if(i != j && triangles[i].vertices[v].x == triangles[j].vertices[u].x && triangles[i].vertices[v].y == triangles[j].vertices[u].y && triangles[i].vertices[v].z == triangles[j].vertices[u].z) {
						if (acos(dot(normalize(triangles[i].normal), normalize(triangles[j].normal))/(length(triangles[i].normal)*length(triangles[j].normal))) < pi/4) {
							vertex = vertex + triangles[j].normal;
							count = count + 1;
						}
					}
				}
			}
			vertex = vertex / float(count);
			triangles[i].normals[v] = normalize(vertex);
		}
	}
	// return triangles;
}

vector<ModelTriangle> parse_obj(string filename, float scale, unordered_map<string, Colour> colours) {

	vector<ModelTriangle> triangles;
	vector<vec3> vertices;
	vector<vec3> normals;
	vector<TexturePoint> texture_points;
	string colour;
	string texture_name;
	bool mirror = false;
	bool refract = false;

	cout << "Texture points: " << texture_points.size() << endl;

	ifstream File(filename);
	string line;

	for (auto it = colours.begin(); it != colours.end(); ++it) {
		cout << it->second.name << endl;
	}
	
	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');
		if(tokens[0] == "vn") {
			vec3 normal(stof(tokens[1]), stof(tokens[2]), stof(tokens[3]));
			normals.push_back(normal);
		} else if(tokens[0] == "v") {
			vec3 vertex(stof(tokens[1])*scale, stof(tokens[2])*scale, stof(tokens[3])*scale);
			vertices.push_back(vertex);
		} else if(tokens[0] == "vt") {
			texture_points.push_back(TexturePoint(stof(tokens[1]), stof(tokens[2])));
		} else if(tokens[0] == "f") {
			vector<string> l1 = split(tokens[1], '/');
			vector<string> l2 = split(tokens[2], '/');
			vector<string> l3 = split(tokens[3], '/');
			ModelTriangle triangle(
				vertices[stoi(l1[0])-1], 
				vertices[stoi(l2[0])-1], 
				vertices[stoi(l3[0])-1], 
				colours[colour]);//Colour("texture_logo.ppm", 255, 255, 255));//
			if(!normals.empty()) {
				triangle.normals[0] = normals[stoi(l1[2])-1]; 
				triangle.normals[1] = normals[stoi(l2[2])-1]; 
				triangle.normals[2] = normals[stoi(l3[2])-1];
			}
			triangle.normal = cross(vec3(triangle.vertices[1]-triangle.vertices[0]),vec3(triangle.vertices[2]-triangle.vertices[0]));
			triangle.mirror = mirror;
			triangle.refract = refract;
			if(!texture_points.empty() && l1[1] != "") {
				triangle.texturePoints[0] = texture_points[stoi(l1[1])-1];
				triangle.texturePoints[1] = texture_points[stoi(l2[1])-1];
				triangle.texturePoints[2] = texture_points[stoi(l3[1])-1];
			} 
			triangles.push_back(triangle);
		}  else if(tokens[0] == "usemtl") {
			if(tokens[1] == "Mirror") mirror = true;
			else mirror = false;
			if(tokens[1] == "Glass") {
				refract = true;
				cout << tokens[1] << endl;
			}
			else refract = false;
			colour = tokens[1];
		}
	}
	if(normals.empty()) {
		vertex_normals(triangles);
		cout << "there are no vertex normals in my obj" << endl;
	} else {
		cout << "there WERE vertex normals in my obj" << endl;
	}
	File.close();
	cout << "Texture points: " << texture_points.size() << endl;

	return triangles;
}

unordered_map<string, Colour> parse_mtl(string filename, unordered_map<string, TextureMap> &textures) {
	unordered_map<string, Colour> colours;
	string colour_name;

	ifstream File(filename);
	string line;

	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');
		if(tokens[0] == "newmtl") {
			colour_name = tokens[1];
		} else if(tokens[0] == "Kd") {
			Colour colour(int(stof(tokens[1])*255),int(stof(tokens[2])*255),int(stof(tokens[3])*255));
			colours.insert({colour_name, colour});
		} else if(tokens[0] == "map_Kd") {
			Colour colour = colours[colour_name];
			colour.name = tokens[1];
			textures.insert({tokens[1], TextureMap(tokens[1])});
			colours[colour_name] = colour;
		}
	}
	File.close();
	return colours;
}

void look_at() {
	vec3 forward = normalize(cam - o);
	vec3 right = normalize(cross(vec3(0.0,1.0,0.0), forward));
	vec3 up = normalize(cross(forward, right));

	cam_orientation[0] = right;
	cam_orientation[1] = up;
	cam_orientation[2] = forward;
}

void reset_camera() {
	cam = vec3(0.0,0.0,4.0);
	// light = vec3(0.0, 1.0, 2.0);
	cam_orientation = mat3(vec3(1.0,0.0,0.0),vec3(0.0,1.0,0.0),vec3(0.0,0.0,1.0));
}

mat3 rotation_y(float t) {
	return mat3(vec3( cos(t), 0.0, sin(t)),vec3(    0.0, 1.0,    0.0),vec3(-sin(t), 0.0, cos(t)));
}
mat3 rotation_x(float t) {
	return mat3(vec3( 1.0,    0.0,    0.0),vec3( 0.0, cos(t),-sin(t)),vec3( 0.0, sin(t), cos(t)));
}
mat3 rotation_z(float t) {
	return mat3(vec3( cos(t),-sin(t), 0.0),vec3( sin(t), cos(t), 0.0),vec3(    0.0,    0.0, 1.0));
}

void orbit(bool orb) {
	if(orb) {
		cam = cam * rotation_y(-pi/180);
		look_at();
	}
}

void move_light(float scale, string dim) {
	for(int i = 0; i < lights.size(); i++) {
		if(dim == "x") {
			lights[i].x += scale;
		} else if(dim == "y") {
			lights[i].y += scale;
		} else if(dim == "z") {
			lights[i].z += scale;
		} else cout << "dimension not known" << endl;
	}
}

function<void(vector<ModelTriangle> &, DrawingWindow &, TextureMap &)> drawing = draw_raytrace;

void animate(vector<ModelTriangle> &triangles, DrawingWindow &window, TextureMap &textures) {
	drawing = draw_wireframe;
	// drawing = draw_rasterise;
	// drawing = draw_raytrace;

	vector<float> steps_24 = {
		3.5000e-01, 2.2208e-01, 1.4091e-01, 8.9405e-02, 5.6728e-02, 3.5994e-02,
        2.2838e-02, 1.4491e-02, 9.1944e-03, 5.8338e-03, 3.7016e-03, 2.3486e-03,
        1.4902e-03, 9.4554e-04, 5.9995e-04, 3.8067e-04, 2.4153e-04, 1.5325e-04,
        9.7239e-05, 6.1698e-05, 3.9148e-05, 2.4839e-05, 1.5760e-05, 1.0000e-05};
	vector<float> steps_48 = {
		2.5000e-02, 2.1166e-02, 1.7920e-02, 1.5172e-02, 1.2846e-02, 1.0876e-02,
        9.2079e-03, 7.7959e-03, 6.6004e-03, 5.5882e-03, 4.7312e-03, 4.0057e-03,
        3.3914e-03, 2.8713e-03, 2.4310e-03, 2.0582e-03, 1.7426e-03, 1.4754e-03,
        1.2491e-03, 1.0576e-03, 8.9538e-04, 7.5807e-04, 6.4182e-04, 5.4340e-04,
        4.6007e-04, 3.8952e-04, 3.2978e-04, 2.7921e-04, 2.3639e-04, 2.0014e-04,
        1.6945e-04, 1.4346e-04, 1.2146e-04, 1.0284e-04, 8.7067e-05, 7.3716e-05,
        6.2411e-05, 5.2840e-05, 4.4737e-05, 3.7877e-05, 3.2068e-05, 2.7151e-05,
        2.2987e-05, 1.9462e-05, 1.6477e-05, 1.3951e-05, 1.1811e-05, 1.0000e-05};

	int n_zero = 5;
	int frames = 0;
	for(int i = 0; i < 48; i++) {
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;
		if(i < 6) {
			cam.x -= 0.05;
		}
		else if(i < 18) {
			cam.x += 0.05;
		}
		else if(i < 24) {
			cam.x -= 0.05;
		}
		frames++;
	}
	drawing = draw_rasterise;
	for(int i = 0; i < 48; i++) {
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;
		if(i < 6) {
			cam.x -= 0.1;
		}
		else if(i < 18) {
			cam.x += 0.1;
		}
		else if(i < 24) {
			cam.x -= 0.1;
		}
		frames++;
	}
	for(int i = 0; i < 48; i++) {
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;
		if(i < 12) {
			cam.z -= 0.1;
		}
		else if(i < 24) {
			cam.z += 0.1;
		}
		else if(i < 36) {
			cam.y += 0.2;
		}
		else if(i < 48) {
			cam.y -= 0.2;
		}
		look_at();
		frames++;
	}
	drawing = draw_raytrace;
	for(int i = 0; i < 96; i++) {
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;
		if(i < 24) {
			move_light(0.075,"x");
			move_light(-0.1,"z");
		}
		else if(i < 48) {
			move_light(-0.15,"x");
		}
		else if(i < 72) {
			move_light(0.075,"x");
			move_light(0.1,"z");
			move_light(-0.075,"y");
		}
		else if(i < 96) {
			move_light(0.075,"y");
		}
		frames++;
	}
	drawing = draw_wireframe;
	for(int i = 0; i < 96; i++) { //1
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		move_light(-0.01,"z");
		float temp_c = i+1;
		float temp_c1 = 96-i;
		// cout << temp_c << endl;
		cam.z -= 0.25*(1/temp_c1);//0.0125;
		cam.x += 0.25*(1/temp_c);//steps_48[i];//temp_c;
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
	for(int i = 0; i < 96; i++) { //2
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		move_light(-0.01,"z");
		float temp_c = 96-i;
		float temp_c1 = i+1;
		// cout << temp_c << endl;
		cam.z -= 0.25*(1/temp_c1);//0.0125;
		cam.x -= 0.25*(1/temp_c);
		cam.y += 0.00625;
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
	for(int i = 0; i < 96; i++) { //3
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		move_light(-0.01,"z");
		float temp_c = i+1;
		float temp_c1 = 96-i;
		// cout << temp_c << endl;
		cam.z -= 0.25*(1/temp_c1);//0.0125;
		cam.x -= 0.25*(1/temp_c);
		cam.y += 0.00625;
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
	for(int i = 0; i < 96; i++) { //4
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		move_light(0.01,"z");
		float temp_c = 96-i;
		float temp_c1 = i+1;
		// cout << temp_c << endl;
		cam.z -= 0.25*(1/temp_c1);
		cam.x += 0.25*(1/temp_c);//0.35*(1/temp_c);
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
	// drawing = draw_rasterise;
	for(int i = 0; i < 96; i++) { //5
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		move_light(0.01,"z");
		// cam.z -= 0.05;
		float temp_c = i+1;
		float temp_c1 = 96-i;
		// cout << temp_c << endl;
		cam.z += 0.25*(1/temp_c1);
		cam.x += 0.25*(1/temp_c);//0.35*(1/temp_c);
		cam.y -= 0.00625;
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
	// drawing = draw_wireframe;
	for(int i = 0; i < 96; i++) {
		draw(window);
		drawing(triangles, window, textures);
		string name = string(n_zero - to_string(frames).length(), '0') + to_string(frames);
		window.savePPM("output/"+name+".ppm");
		cout << "saved " << frames << endl;

		cam.y -= 0.00625;
		float temp_c = 96-i; 
		// float temp_c = i+1;
		// temp_c = 0.35*(1/temp_c);
		cam.z += 0.25*(1/temp_c);;
		// cout << temp_c << endl;
		// cam.x += temp_c;
		// cam = cam * rotation_y(pi/180);
		look_at();
		frames++;
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_PAGEDOWN) cam.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_PAGEUP) cam.y += 0.1;
		else if (event.key.keysym.sym == SDLK_w) {cam.z -= 0.1; cout << "[" << cam.x << "," << cam.y << "," << cam.z << "]" << endl;}
		else if (event.key.keysym.sym == SDLK_a) cam.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_s) cam.z += 0.1;
		else if (event.key.keysym.sym == SDLK_d) cam.x += 0.1;
		else if (event.key.keysym.sym == SDLK_q)     cam = cam * rotation_y(-pi/180);
		else if (event.key.keysym.sym == SDLK_e)     cam = cam * rotation_y( pi/180);
		else if (event.key.keysym.sym == SDLK_EQUALS)     cam = cam * rotation_x(-pi/180);
		else if (event.key.keysym.sym == SDLK_MINUS)     cam = cam * rotation_x( pi/180);
		else if (event.key.keysym.sym == SDLK_LEFT)  cam_orientation = cam_orientation * rotation_y(-pi/180);
		else if (event.key.keysym.sym == SDLK_RIGHT) cam_orientation = cam_orientation * rotation_y( pi/180);
		else if (event.key.keysym.sym == SDLK_UP)    cam_orientation = cam_orientation * rotation_x(-pi/180);
		else if (event.key.keysym.sym == SDLK_DOWN)  cam_orientation = cam_orientation * rotation_x( pi/180);
		else if (event.key.keysym.sym == SDLK_z)     cam = cam * rotation_z(-pi/180);
		else if (event.key.keysym.sym == SDLK_x)     cam = cam * rotation_z( pi/180);
		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		else if (event.key.keysym.sym == SDLK_l) look_at();
		else if (event.key.keysym.sym == SDLK_r) reset_camera();
		else if (event.key.keysym.sym == SDLK_1) { drawing = draw_raytrace; cout << "[drawing]: raytrace" << endl; }
		else if (event.key.keysym.sym == SDLK_2) { drawing = draw_rasterise; cout << "[drawing]: rasterise" << endl; }
		else if (event.key.keysym.sym == SDLK_3) { drawing = draw_wireframe; cout << "[drawing]: wireframe" << endl; }
		else if (event.key.keysym.sym == SDLK_4) { brightness = get_scale; cout << "[lighting]: scale" << endl; }
		else if (event.key.keysym.sym == SDLK_5) { brightness = gourad; cout << "[lighting]: gourad" << endl; }
		else if (event.key.keysym.sym == SDLK_6) { brightness = phong; cout << "[lighting]: phong" << endl; }
		else if (event.key.keysym.sym == SDLK_KP_8) move_light(-0.1, "z");
		else if (event.key.keysym.sym == SDLK_KP_2) move_light(0.1, "z");
		else if (event.key.keysym.sym == SDLK_KP_6) move_light(0.1, "x");
		else if (event.key.keysym.sym == SDLK_KP_4) move_light(-0.1, "x");
		else if (event.key.keysym.sym == SDLK_KP_MINUS) move_light(-0.1, "y");
		else if (event.key.keysym.sym == SDLK_KP_PLUS) move_light(0.1, "y");
		else if (event.key.keysym.sym == SDLK_LEFTBRACKET)  { proximity = (proximity) ? false : true; cout << "[proximity]: " << proximity << endl; }
		else if (event.key.keysym.sym == SDLK_RIGHTBRACKET) { angle_of  = (angle_of)  ? false : true; cout << "[angle_of]: " << angle_of << endl; }
		else if (event.key.keysym.sym == SDLK_HASH)         { shadows   = (shadows)   ? false : true; cout << "[shadows]: " << shadows << endl; }
		else if (event.key.keysym.sym == SDLK_QUOTE)        { specular  = (specular)  ? false : true; cout << "[specular]: " << specular << endl; }
		else if (event.key.keysym.sym == SDLK_p) show_light = (show_light) ? false : true;
		else if (event.key.keysym.sym == SDLK_7) {window.savePPM("output.ppm"); cout << "[SAVED]" << endl;}
		else if (event.key.keysym.sym == SDLK_7) {; cout << "[SAVED]" << endl;}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {window.savePPM("output.ppm"); cout << "[SAVED]" << endl;}
}

int main(int argc, char *argv[]) {

	time_t start,finish;
	unordered_map<string, TextureMap> textures;
	// vector<ModelTriangle> t_0 = parse_obj("logo2.obj", 0.002, parse_mtl("logo.mtl"));
	// vector<ModelTriangle> t = parse_obj("low_poly_bunny.obj", 0.5, parse_mtl("cornell-box.mtl", textures));
	// vector<ModelTriangle> t = parse_obj("sphere.obj", 0.5, parse_mtl("cornell-box.mtl", textures));
	vector<ModelTriangle> t = parse_obj("cornell-box.obj", 0.5, parse_mtl("cornell-box.mtl", textures));
	TextureMap texture("logo.ppm");
	// vector<ModelTriangle> t = parse_obj("cornell-rabbit.obj", 0.5, parse_mtl("cornell-box.mtl", textures));
	// vector<ModelTriangle> t_2 = parse_obj("bunny.obj", 0.01, parse_mtl("cornell-box.mtl"));
	// for(int i = 0; i < t.size(); i++) {
	// 	cout << t[i] << endl;
	// }
	// t.insert(t.end(), t_2.begin(), t_2.end());
	// t.insert(t.end(), t_0.begin(), t_0.end());
	cout << "Triangles: " << t.size() << endl;

	
	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	animate(t, window_grey, texture);
	draw(window_grey);
	time(&start);
	drawing(t, window_grey, texture);
	time(&finish);
	cout << "rendering took: " << difftime(finish,start) << " seconds" << endl;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		orbit(orbiting);

		draw(window_grey);
		drawing(t, window_grey, texture);
		// window_grey.savePPM("output.ppm"); 
		// cout << "[SAVED]" << endl;
		// drawing(t_2, window_grey);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}
