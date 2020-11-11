#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp> 
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

#define WIDTH 600
#define HEIGHT 600

#define PI 3.14159265359

glm::vec3 camera(0.0, 0.0, 4.0);
float distance = 700;
glm::mat3 cameraOrientation(
	glm::vec3(1.0, 0.0, 0.0),
	glm::vec3(0.0, 1.0, 0.0),
	glm::vec3(0.0, 0.0, 1.0)
);

bool orbiting = false;

std::vector<float> interpolateSingleFloats(float from, float to, int numVals) {
	std::vector<float> result;
	float step = (to - from)/(numVals-1); 
	float temp = from; 

	result.push_back(temp);

	for (int i = 0; i < numVals-1; i++) {
		temp = temp + step;
		result.push_back(temp);
	}
	return result;
}

std::vector<CanvasPoint> interpolatePoints(CanvasPoint start, CanvasPoint end, int steps){
	std::vector<CanvasPoint> result;
	float stepX = (end.x - start.x)/(steps-1);
	float stepY = (end.y - start.y)/(steps-1);
	float stepDepth = (end.depth - start.depth)/(steps-1);
	
	CanvasPoint temp = start;
	result.push_back(temp);

	for (int i = 0; i < steps-1; i++) {
		temp.x = temp.x + stepX;
		temp.y = temp.y + stepY;
		temp.depth = temp.depth + stepDepth;

		result.push_back(temp);
	}
	return result;
}

std::vector<TexturePoint> interpolatePoints(TexturePoint start, TexturePoint end, int steps){
	std::vector<TexturePoint> result;
	float stepX = (end.x - start.x)/(steps-1);
	float stepY = (end.y - start.y)/(steps-1);

	TexturePoint temp = start;
	result.push_back(temp);

	for (int i = 0; i < steps-1; i++) {
		temp.x = temp.x + stepX;
		temp.y = temp.y + stepY;
		result.push_back(temp);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numVals) {
	std::vector<glm::vec3> result;
	float xStep = (to.x - from.x)/(numVals - 1);
	float yStep = (to.y - from.y)/(numVals - 1);
	float zStep = (to.z - from.z)/(numVals - 1);

	float xTemp = from.x;
	float yTemp = from.y;
	float zTemp = from.z;

	glm::vec3 first(xTemp, yTemp, zTemp);
	result.push_back(first);

	for (int i = 0; i < numVals-1; i++) {
		xTemp = xTemp + xStep;
		yTemp = yTemp + yStep;
		zTemp = zTemp + zStep;
		glm::vec3 temp(xTemp, yTemp, zTemp); 
		result.push_back(temp);
	}

	return result;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numSteps;
	float yStepSize = yDiff/numSteps;
	for (float i = 0.0; i < numSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
		window.setPixelColour(round(x), round(y), set);
	}
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void textureFill(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.texturePoint, mid.texturePoint);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
		std::swap(mid.texturePoint, top.texturePoint);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.texturePoint, mid.texturePoint);
	}
	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	// TOP TRIANGLE ---------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<CanvasPoint> left = interpolatePoints(top, mid, mid.y-top.y+2);
	std::vector<TexturePoint> leftTexture = interpolatePoints(top.texturePoint, mid.texturePoint, mid.y-top.y+2);

	std::vector<CanvasPoint> right = interpolatePoints(top, split, mid.y-top.y+2);
	std::vector<TexturePoint> rightTexture = interpolatePoints(top.texturePoint, split.texturePoint, mid.y-top.y+2);

	for (int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);

		std::vector<TexturePoint> texturePoints = interpolatePoints(leftTexture[i], rightTexture[i], steps+2);

		for (int c = 0; c < texturePoints.size(); c++) {
			int x_coord = texturePoints.at(c).x;
			int y_coord = texturePoints.at(c).y;
			uint32_t col = texture.pixels.at((y_coord*texture.width) + x_coord);

			window.setPixelColour(round(points[c].x), round(points[i].y), col);
		}

	}

	// BOTTOM TRIANGLE ------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<CanvasPoint> left2 = interpolatePoints(bot, mid, bot.y-mid.y+2);
	std::vector<TexturePoint> leftTexture2 = interpolatePoints(bot.texturePoint, mid.texturePoint, bot.y-mid.y+2);

	std::vector<CanvasPoint> right2 = interpolatePoints(bot, split, bot.y-mid.y+2);
	std::vector<TexturePoint> rightTexture2 = interpolatePoints(bot.texturePoint, split.texturePoint, bot.y-mid.y+2);

	for (int i = 0; i < left2.size(); i++) {
		int steps = abs(left2[i].x - right2[i].x);
				
		std::vector<CanvasPoint> points = interpolatePoints(left2[i], right2[i], steps+2);

		std::vector<TexturePoint> texturePoints = interpolatePoints(leftTexture2[i], rightTexture2[i], steps+2);

		for (int c = 0; c < texturePoints.size(); c++) {
			int x_coord = texturePoints.at(c).x;
			int y_coord = texturePoints.at(c).y;
			uint32_t col = texture.pixels.at(int((y_coord*texture.width) + x_coord));

			window.setPixelColour(round(points[c].x), round(points[i].y), col);
		}

	}

}

void fillCornell(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.depth, mid.depth);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
		std::swap(mid.depth, top.depth);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.depth, mid.depth);
	}
	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	split.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	// TOP TRIANGLE ---------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<CanvasPoint> left = interpolatePoints(top, mid, mid.y-top.y+2);

	std::vector<CanvasPoint> right = interpolatePoints(top, split, mid.y-top.y+2);

	for (int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);
		
		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height){
				if (-1/points[c].depth > depths[newX][newY]) {

					depths[newX][newY] = -1/points[c].depth;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(newX, newY, set);
				}
			}
			

		}

	}

	// BOTTOM TRIANGLE ------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<CanvasPoint> left2 = interpolatePoints(bot, mid, bot.y-mid.y+2);

	std::vector<CanvasPoint> right2 = interpolatePoints(bot, split, bot.y-mid.y+2);

	for (int i = 0; i < left2.size(); i++) {
		int steps = abs(left2[i].x - right2[i].x);
				
		std::vector<CanvasPoint> points = interpolatePoints(left2[i], right2[i], steps+2);

		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height){
				if (-1/points[c].depth > depths[newX][newY]) {

					depths[newX][newY] = -1/points[c].depth;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(newX, newY, set);
				}
			}
			

		}

	}

}

void drawCornellWireframe(DrawingWindow &window, std::vector<ModelTriangle> triangles) {
	for (int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		for (int j = 0; j < 3; j++) {
			int u = -(distance * triangles[i].vertices[j].x/(triangles[i].vertices[j].z - camera.z)) + (window.width / 2);
			int v = (distance * triangles[i].vertices[j].y/(triangles[i].vertices[j].z - camera.z)) + (window.height / 2);

			triangle.vertices[j] = CanvasPoint(u, v);
		}
		
		drawTriangle(window, triangle, Colour(255,255,255));
 	}

}



void drawCornell(DrawingWindow &window, std::vector<ModelTriangle> triangles) {
	std::vector<std::vector<float>> depths(window.width, std::vector<float> (window.height, -std::numeric_limits<float>::infinity()));

	for (int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		bool isTexture = false;
		TextureMap texture;

		if (triangles[i].colour.name != "") {
			texture = TextureMap(triangles[i].colour.name);
			isTexture = true;
		}
		for (int j = 0; j < 3; j++) {
			glm::vec3 cameraToVertex = glm::vec3(triangles[i].vertices[j].x - camera.x, triangles[i].vertices[j].y - camera.y, triangles[i].vertices[j].z - camera.z);

			glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;

			int u = -(distance * (adjustedVector.x)/(adjustedVector.z)) + (window.width / 2);
			int v = (distance * (adjustedVector.y)/(adjustedVector.z)) + (window.height / 2);

			triangle.vertices[j] = CanvasPoint(u, v, adjustedVector.z);

			if (isTexture == true) {
				triangle.vertices[j].texturePoint = triangles[i].texturePoints[j];
				triangle.vertices[j].texturePoint.x *= texture.width;
				triangle.vertices[j].texturePoint.y *= texture.height;

			}	
		}

		//std::cout << triangles[i].colour.name << std::endl;

		if (isTexture == true) {
			textureFill(window, triangle, texture);
		} else fillCornell(window, triangle, triangles[i].colour, depths);

		
		
 	}

}

void lookAt() {
	glm::vec3 forward = glm::normalize(camera - glm::vec3(0,0,0));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0,1,0), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
}

void resetCamera() {
	camera[0] = 0.0;
	camera[1] = 0.0;
	camera[2] = 4.0;

	cameraOrientation[0] = glm::vec3(1.0, 0.0, 0.0);
	cameraOrientation[1] = glm::vec3(0.0, 1.0, 0.0);
	cameraOrientation[2] = glm::vec3(0.0, 0.0, 1.0);
}

void orbit(bool orb) {
	if(orb) {
		float theta = -PI/180;
			glm::mat3 m = glm::mat3(
				cos(theta), 0, sin(theta),
				0, 1, 0,
				-sin(theta), 0, cos(theta)
			);
			camera = camera * m;
			lookAt();
	}
}

std::vector<ModelTriangle> parseObj(std::string filename, float scale, std::unordered_map<std::string, Colour> colours) {
	std::vector<ModelTriangle> output;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> textureVertices;
	std::string colour;

	std::ifstream File(filename);
	std::string line;
	
	while(std::getline(File, line)) {
		if(line == "") continue;

		std::vector<std::string> tokens = split(line, ' ');

		if (tokens[0] == "v") {
			glm::vec3 temp = glm::vec3(stof(tokens[1])*scale, stof(tokens[2])*scale, stof(tokens[3])*scale); 
			vertices.push_back(temp);

		} else if (tokens[0] == "f") {
			// when no texture map vertices, the second vector item is equal to ""
			std::vector<std::string> a = split(tokens[1],'/');
			std::vector<std::string> b = split(tokens[2],'/');
			std::vector<std::string> c = split(tokens[3],'/');

			if (a[1] == "") {
				ModelTriangle triangle(vertices[stoi(a[0])-1], vertices[stoi(b[0])-1], vertices[stoi(c[0])-1], colours[colour]);
				output.push_back(triangle);
			} else {
				ModelTriangle triangle(vertices[stoi(a[0])-1], vertices[stoi(b[0])-1], vertices[stoi(c[0])-1], colours[colour]);
				triangle.texturePoints[0] = textureVertices[stoi(a[1])-1];
				triangle.texturePoints[1] = textureVertices[stoi(b[1])-1];
				triangle.texturePoints[2] = textureVertices[stoi(c[1])-1];
				output.push_back(triangle);
			}
		} else if (tokens[0] == "usemtl") {
			colour = tokens[1];
		} else if (tokens[0] == "vt") {
			TexturePoint temp = TexturePoint(stof(tokens[1]), stof(tokens[2]));
			textureVertices.push_back(temp);
		}
	}

	File.close();

	return output;
}

std::unordered_map<std::string, Colour> parseMtl(std::string filename) {
	std::unordered_map<std::string, Colour> colours;
	std::string colour;

	std::ifstream File(filename);
	std::string line;

	while(std::getline(File, line)) {
		if(line == "") continue;

		std::vector<std::string> tokens = split(line, ' ');

		if (tokens[0] == "newmtl") {
			colour = tokens[1];
		} else if (tokens[0] == "Kd") {
			std::string a = tokens[1];
			std::string b = tokens[2];
			std::string c = tokens[3];

			Colour temp(int(stof(a)*255), int(stof(b)*255), int(stof(c)*255));
			colours.insert({colour, temp});
		} else if (tokens[0] == "map_Kd") {
			Colour temp = colours[colour];
			temp.name = tokens[1];
			colours[colour] = temp;
		}
	}

	return colours;
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
		}
	}
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) camera.x += 0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT) camera.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_UP) camera.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_DOWN) camera.y += 0.1;
		else if (event.key.keysym.sym == SDLK_s) camera.z += 0.1;
		else if (event.key.keysym.sym == SDLK_w) camera.z -= 0.1;
		// CAMERA ROTATION
		else if (event.key.keysym.sym == SDLK_r) {
			float theta = -PI/180;
			glm::mat3 m = glm::mat3(
				1, 0, 0,
				0, cos(theta), -sin(theta),
				0, sin(theta), cos(theta)
			);
			camera = camera * m;
			lookAt();
		}
		else if (event.key.keysym.sym == SDLK_e) {
			float theta = PI/180;
			glm::mat3 m = glm::mat3(
				1, 0, 0,
				0, cos(theta), -sin(theta),
				0, sin(theta), cos(theta)
			);
			camera = camera * m;
			lookAt();
		}
		else if (event.key.keysym.sym == SDLK_f) {
			float theta = PI/180;
			glm::mat3 m = glm::mat3(
				cos(theta), 0, sin(theta),
				0, 1, 0,
				-sin(theta), 0, cos(theta)
			);
			camera = camera * m;
			lookAt();
		}
		else if (event.key.keysym.sym == SDLK_g) {
			float theta = -PI/180;
			glm::mat3 m = glm::mat3(
				cos(theta), 0, sin(theta),
				0, 1, 0,
				-sin(theta), 0, cos(theta)
			);
			camera = camera * m;
			lookAt();
		}
		// CAMERA ORIENTATION ROTATION
		else if (event.key.keysym.sym == SDLK_u) {
			float theta = PI/180;
			glm::mat3 m = glm::mat3(
				glm::vec3(1, 0, 0),
				glm::vec3(0, cos(theta), -sin(theta)),
				glm::vec3(0, sin(theta), cos(theta))
			);
			cameraOrientation = cameraOrientation * m;
		}
		else if (event.key.keysym.sym == SDLK_j) {
			float theta = -PI/180;
			glm::mat3 m = glm::mat3(
				glm::vec3(1, 0, 0),
				glm::vec3(0, cos(theta), -sin(theta)),
				glm::vec3(0, sin(theta), cos(theta))
			);
			cameraOrientation = cameraOrientation * m;
		}
		else if (event.key.keysym.sym == SDLK_k) {
			float theta = -PI/180;
			glm::mat3 m = glm::mat3(
				glm::vec3(cos(theta), 0, sin(theta)),
				glm::vec3(0, 1, 0),
				glm::vec3(-sin(theta), 0, cos(theta))
			);
			cameraOrientation = cameraOrientation * m;
		}
		else if (event.key.keysym.sym == SDLK_h) {
			float theta = PI/180;
			glm::mat3 m = glm::mat3(
				glm::vec3(cos(theta), 0, sin(theta)),
				glm::vec3(0, 1, 0),
				glm::vec3(-sin(theta), 0, cos(theta))
			);
			cameraOrientation = cameraOrientation * m;
		}
		else if (event.key.keysym.sym == SDLK_q) {
			lookAt();
		}
		else if (event.key.keysym.sym == SDLK_1) {
			resetCamera();
		}
		else if (event.key.keysym.sym == SDLK_o) {
			orbiting = (orbiting) ? false : true;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	std::vector<ModelTriangle> triangles;
	std::unordered_map<std::string, Colour> colours;

	colours = parseMtl("cornell-box.mtl");
	triangles = parseObj("cornell-box.obj", 0.4, colours);

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		//update(window);
		orbit(orbiting);
		draw(window);
		drawCornell(window, triangles);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
