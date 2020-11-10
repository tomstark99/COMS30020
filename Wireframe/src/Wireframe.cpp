#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

using namespace std;
using namespace glm;

#define WIDTH 600
#define HEIGHT 600

vec3 cam(0.0, 0.0, 4.0);
float focal = 500.0;

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

std::vector<CanvasPoint> interpolate_canvas(CanvasPoint from, CanvasPoint to, int steps){
	std::vector<CanvasPoint> points;
	float x_step = (from.x - to.x)/(steps-1);
	float y_step = (from.y - to.y)/(steps-1);

	CanvasPoint temp = from;
	points.push_back(temp);

	for (int i = 0; i < steps-1; i++) {
		temp.x = temp.x + x_step;
		temp.y = temp.y + y_step;
		points.push_back(temp);
	}
	return points;
}

std::vector<TexturePoint> interpolate_texture(TexturePoint from, TexturePoint to, int steps){
	std::vector<TexturePoint> points;
	float x_step = (from.x - to.x)/(steps-1);
	float y_step = (from.y - to.y)/(steps-1);

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

void draw_line_depths(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window, vector<vector<float>> &depths) {
	float x_diff = to.x - from.x;
	float y_diff = to.y - from.y;
	float d_diff = to.depth - from.depth;
	float steps = std::max(abs(x_diff),abs(y_diff));

	float x_step = x_diff/steps;
	float y_step = y_diff/steps;
	float d_step = d_diff/steps;

	uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

	for(float i = 0.0; i < steps; i++) {
		float x = from.x + (x_step * i);
		float y = from.y + (y_step * i);
		float d = from.depth + (d_step * i);
		if(-1/d > depths[round(x)][round(y)]) {
			depths[round(x)][round(y)] = -1/d;
			window.setPixelColour(round(x), round(y), c);
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

	float x1_diff = mid.x - top.x;
	float y1_diff = mid.y - top.y;
	float x2_diff = bot.x - top.x;
	float y2_diff = bot.y - top.y;
	float d1_diff = mid.depth - top.depth;
	float d2_diff = bot.depth - top.depth;

	float steps1 = std::max(abs(x1_diff),abs(y1_diff));
	float steps2 = std::max(abs(x2_diff),abs(y2_diff));

	float x1_step = x1_diff/steps1;
	float x2_step = x2_diff/steps2;
	float y1_step = y1_diff/steps1;
	float y2_step = y2_diff/steps2;
	float d1_step = d1_diff/steps1;
	float d2_step = d2_diff/steps2;

	for(float i = 0.0; i < steps1; i++) {
		for(float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1_step * i);
			float y1 = top.y + (y1_step * i);
			float x2 = top.x + (x2_step * j);
			float y2 = top.y + (y2_step * j);
			float d1 = top.depth + (d1_step * i);
			float d2 = top.depth + (d2_step * j);
			draw_line_depths(CanvasPoint(round(x1),round(y1),d1),CanvasPoint(round(x2),round(y2),d2), colour, window, depths);
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
	CanvasTriangle t_2 = CanvasTriangle(mid_2,mid,bot);

	fill_half_triangle(t_1, colour, window, depths);
	fill_half_triangle(t_2, colour, window, depths);

	// draw_triangle(triangle, Colour(255,255,255), window);
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

void draw_obj(vector<ModelTriangle> triangles, DrawingWindow &window) {

	vector<vector<float>> depths(window.width, vector<float>(window.height, -(numeric_limits<float>::infinity())));

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		CanvasTriangle t;
		for(int i = 0; i < triangle.vertices.size(); i++) {
			vec3 vertex = triangle.vertices[i];
			int u = -(focal * vertex.x)/(vertex.z-cam.z) + window.width/2;
			int v = (focal * vertex.y)/(vertex.z-cam.z) + window.height/2;
			t.vertices[i] = CanvasPoint(u,v, (vertex.z-cam.z));
		}
		// draw_triangle(t, Colour(255,255,255), window);
		fill_triangle(t, triangle.colour, window, depths);
	}
}

vector<ModelTriangle> parse_obj(string filename, float scale, unordered_map<string, Colour> colours) {

	vector<ModelTriangle> triangles;
	vector<vec3> vertices;
	string colour;

	ifstream File(filename);
	string line;
	
	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');
		if(tokens[0] == "v") {
			vec3 vertex(stof(tokens[1])*scale, stof(tokens[2])*scale, stof(tokens[3])*scale);
			vertices.push_back(vertex);
		} else if(tokens[0] == "f") {
			ModelTriangle triangle(vertices[stoi(tokens[1])-1], vertices[stoi(tokens[2])-1], vertices[stoi(tokens[3])-1], colours[colour]);
			triangles.push_back(triangle);
		} else if(tokens[0] == "usemtl") {
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
		}
	}
	File.close();
	return colours;
}

int main(int argc, char *argv[]) {

	vector<ModelTriangle> t = parse_obj("cornell-box.obj", 0.5, parse_mtl("cornell-box.mtl"));

	// for(int i = 0; i < t.size(); i++) {
	// 	cout << t[i].colour.red << ' ' << t[i].colour.green << ' ' << t[i].colour.blue << endl;
	// }

	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		update(window_grey);
		draw_obj(t, window_grey);
		// texture_triangle(texture, t, window_grey);
		// fill_triangle(t, Colour(0,0,255),window_grey);
		// draw_triangle(t, Colour(255,255,255), window_grey);
		// random_triangle(window_grey);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}
