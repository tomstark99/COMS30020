#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <utility>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

using namespace std;
using namespace glm;

#define WIDTH 600
#define HEIGHT 600

#define pi 3.14159265359

vec3 cam(0.0, 0.0, 4.0);
vec3   o(0.0, 0.0, 0.0);
mat3 cam_orientation(
	vec3(1.0,0.0,0.0),
	vec3(0.0,1.0,0.0),
	vec3(0.0,0.0,1.0)
);
float focal = 500.0;
bool orbiting = false;

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
		float x = from.x + (x_step * i);
		float y = from.y + (y_step * i);
		window.setPixelColour(round(x), round(y), c);
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
					// uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
					window.setPixelColour(x, y, texture.pixels[round(points_texture[j].y)*texture.width + round(points_texture[j].x)]);
					// window.setPixelColour(x,y,c);
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

void draw_obj(vector<ModelTriangle> triangles, DrawingWindow &window) {

	vector<vector<float>> depths(window.width, vector<float>(window.height, -(numeric_limits<float>::infinity())));

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		CanvasTriangle t;
		for(int j = 0; j < triangle.vertices.size(); j++) {
			vec3 vertex = triangle.vertices[j];
			vec3 cam_to_vertex(vertex.x - cam.x, vertex.y - cam.y, vertex.z - cam.z);
			vec3 adjusted_vertex = cam_to_vertex * cam_orientation;
		
			int u = -(focal * (adjusted_vertex.x)/(adjusted_vertex.z)) + (window.width/2);
			int v = (focal * (adjusted_vertex.y)/(adjusted_vertex.z)) + (window.height/2);
			;
			t.vertices[j] = CanvasPoint(u,v, adjusted_vertex.z);
			t.vertices[j].texturePoint = triangle.texturePoints[j];
		}
		// draw_triangle(t, Colour(255,255,255), window);
		if (triangle.colour.name != "") {
			TextureMap texture(triangle.colour.name);
			for(int j = 0; j < t.vertices.size(); j++) {
				t.vertices[j].texturePoint.x *= texture.width;
				t.vertices[j].texturePoint.y *= texture.height;
			}
			texture_triangle(texture, t, window, depths);
		} else {
			fill_triangle(t, triangle.colour, window, depths);
		}
	}
}

vector<ModelTriangle> parse_obj(string filename, float scale, unordered_map<string, Colour> colours) {

	vector<ModelTriangle> triangles;
	vector<vec3> vertices;
	vector<TexturePoint> texture_points;
	string colour;
	string texture_name;

	ifstream File(filename);
	string line;
	
	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');
		if(tokens[0] == "v") {
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
				colours[colour]);
			if(l1[1] != "") {
				triangle.texturePoints[0] = texture_points[stoi(l1[1])-1];
				triangle.texturePoints[1] = texture_points[stoi(l2[1])-1];
				triangle.texturePoints[2] = texture_points[stoi(l3[1])-1];
			} 
			triangles.push_back(triangle);
		}  else if(tokens[0] == "usemtl") {
			if(tokens[1] != "") {
				texture_name = tokens[1];
			}
			colour = tokens[1];
		}
	}
	File.close();
	return triangles;
}

unordered_map<string, Colour> parse_mtl(string filename) {
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
			colours[colour_name] = colour;
		}
	}
	File.close();
	return colours;
}

void look_at() {
	vec3 forward = normalize(cam - vec3(0.0,0.0,0.0));
	vec3 right = normalize(cross(vec3(0.0,1.0,0.0), forward));
	vec3 up = normalize(cross(forward, right));

	cam_orientation[0] = right;
	cam_orientation[1] = up;
	cam_orientation[2] = forward;
}

void reset_camera() {
	cam = vec3(0.0,0.0,4.0);
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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_PAGEDOWN) cam.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_PAGEUP) cam.y += 0.1;
		else if (event.key.keysym.sym == SDLK_w) cam.z -= 0.1;
		else if (event.key.keysym.sym == SDLK_a) cam.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_s) cam.z += 0.1;
		else if (event.key.keysym.sym == SDLK_d) cam.x += 0.1;
		else if (event.key.keysym.sym == SDLK_q)     cam = cam * rotation_y(-pi/180);
		else if (event.key.keysym.sym == SDLK_e)     cam = cam * rotation_y( pi/180);
		else if (event.key.keysym.sym == SDLK_0)     cam = cam * rotation_x(-pi/180);
		else if (event.key.keysym.sym == SDLK_9)     cam = cam * rotation_x( pi/180);
		else if (event.key.keysym.sym == SDLK_LEFT)  cam_orientation = cam_orientation * rotation_y(-pi/180);
		else if (event.key.keysym.sym == SDLK_RIGHT) cam_orientation = cam_orientation * rotation_y( pi/180);
		else if (event.key.keysym.sym == SDLK_UP)    cam_orientation = cam_orientation * rotation_x(-pi/180);
		else if (event.key.keysym.sym == SDLK_DOWN)  cam_orientation = cam_orientation * rotation_x( pi/180);
		else if (event.key.keysym.sym == SDLK_z)     cam = cam * rotation_z(-pi/180);
		else if (event.key.keysym.sym == SDLK_x)     cam = cam * rotation_z( pi/180);
		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		else if (event.key.keysym.sym == SDLK_l) look_at();
		else if (event.key.keysym.sym == SDLK_r) reset_camera();
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {

	vector<ModelTriangle> t = parse_obj("cornell-box.obj", 0.5, parse_mtl("cornell-box.mtl"));

	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		orbit(orbiting);
		draw(window_grey);
		draw_obj(t, window_grey);
		// texture_triangle(texture, t, window_grey);
		// fill_triangle(t, Colour(0,0,255),window_grey);
		// draw_triangle(t, Colour(255,255,255), window_grey);
		// random_triangle(window_grey);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}
