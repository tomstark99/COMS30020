#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
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
#define WIREFRAME 1
#define RASTERISE 2
#define RAYTRACE 3

#define PI 3.14159265359

glm::vec3 camera(0.0, 0.0, 4.0);
float distance = 700;
glm::mat3 cameraOrientation(
	glm::vec3(1.0, 0.0, 0.0),
	glm::vec3(0.0, 1.0, 0.0),
	glm::vec3(0.0, 0.0, 1.0)
);
int drawing = 1;
glm::vec3 light(0.0,0.85,0.0);
//glm::vec3 light(-0.64901096, 2.739334, 0.532032);
//glm::vec3 light(-0.64901096, 2.7384973, -0.51796794);
//glm::vec3 light(0.650989, 2.7384973, -0.51796794);
//glm::vec3 light(0.650989, 2.739334, 0.532032);

bool proximity = true, angle_of = true, shadows = true;

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

bool inShadow(std::vector<ModelTriangle> triangles, glm::vec3 intersectionPoint, size_t index) {
	bool shadow = false;
	glm::vec3 shadowRay = light - intersectionPoint;
	float length = glm::length(shadowRay);

	for (int i = 0; i < triangles.size(); i++) {
		if (i != index) {
			ModelTriangle triangle = triangles[i];
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			glm::vec3 SPVector = intersectionPoint - triangle.vertices[0];
			glm::mat3 DEMatrix(-normalize(shadowRay), e0, e1);
			glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector; 
			float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

			if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0) && t > 0.01f && t < length) {
				shadow = true;
				break;
			}
		}

	}

	return shadow;
}

float getBrightness(glm::vec3 intersectionPoint, glm::vec3 normal) {
	// glm::vec3 lightRay = light - intersectionPoint;
	glm::vec3 lightRay = intersectionPoint - light;
	float length = glm::length(lightRay);
	float brightness = (proximity) ? 1/(length*length) : 0;

	// lightRay = glm::normalize(lightRay);
	// normal = glm::normalize(normal);

	float angleOfIncidence = glm::dot(lightRay, normal);

	//std::cout << angleOfIncidence << std::endl;

	if (angleOfIncidence > 0 && angle_of) {
		brightness += angleOfIncidence;
	} 

	if (brightness > 1) {
		brightness = 1;
	}

	return brightness;
}

RayTriangleIntersection getClosestIntersection(std::vector<ModelTriangle> triangles, glm::vec3 rayDirection) {
	RayTriangleIntersection closestIntersection;
	closestIntersection.distanceFromCamera = std::numeric_limits<float>::infinity();

	// possibleSolution returns t,u,v
	// t = distance along the ray from the camera to the intersection point
	// u = the proportion along the triangle's first edge that the intersection point occurs
	// v = the proportion along the triangle's second edge that the intersection point occurs
	for (int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = camera - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector; 
		float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && ((u + v) <= 1.0)) {
			if (closestIntersection.distanceFromCamera > t && t > 0) {
				closestIntersection.distanceFromCamera = t;
				closestIntersection.intersectedTriangle = triangle;
				closestIntersection.triangleIndex = i;
				glm::vec3 intersectionPoint = triangle.vertices[0] + 
					u*(triangle.vertices[1]-triangle.vertices[0]) + 
					v*(triangle.vertices[2]-triangle.vertices[0]);

				closestIntersection.intersectionPoint = intersectionPoint;
			}
		}
	}
	return closestIntersection;

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
		if (round(x) >= 0 && round(x) < window.width && round(y) >= 0 && round(y) < window.height) {
			window.setPixelColour(round(x), round(y), set);
		}

	}
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void textureFill(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture, std::vector<std::vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.texturePoint, mid.texturePoint);
		std::swap(bot.depth, mid.depth);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
		std::swap(mid.texturePoint, top.texturePoint);
		std::swap(mid.depth, top.depth);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
		std::swap(bot.texturePoint, mid.texturePoint);
		std::swap(bot.depth, mid.depth);
	}
	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	split.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	// converting depths to inverse for later interpolation
	top.depth = -1/top.depth;
	mid.depth = -1/mid.depth;
	bot.depth = -1/bot.depth;
	split.depth = -1/split.depth;

	// TOP TRIANGLE ---------------------------------------------------------------------------------------------------------------------------------------------------

	std::vector<CanvasPoint> left = interpolatePoints(top, mid, mid.y-top.y+2);
	std::vector<TexturePoint> leftTexture = interpolatePoints(top.texturePoint, mid.texturePoint, mid.y-top.y+2);

	std::vector<CanvasPoint> right = interpolatePoints(top, split, mid.y-top.y+2);
	std::vector<TexturePoint> rightTexture = interpolatePoints(top.texturePoint, split.texturePoint, mid.y-top.y+2);

	for (int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		std::vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);

		std::vector<TexturePoint> texturePoints = interpolatePoints(leftTexture[i], rightTexture[i], steps+2);

		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height) {
				if (points[c].depth > depths[newX][newY]) {
					depths[newX][newY] = points[c].depth;
					int x_coord = texturePoints.at(c).x;
					int y_coord = texturePoints.at(c).y;
					uint32_t col = texture.pixels.at(int((y_coord*texture.width) + x_coord));

					window.setPixelColour(newX, newY, col);	
				}			
			}
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

		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height) {
				if (points[c].depth > depths[newX][newY]) {
					depths[newX][newY] = points[c].depth;
					int x_coord = texturePoints.at(c).x;
					int y_coord = texturePoints.at(c).y;
					uint32_t col = texture.pixels.at(int((y_coord*texture.width) + x_coord));

					window.setPixelColour(newX, newY, col);	
				}			
			}

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

	// converting depths to inverse for later interpolation
	top.depth = -1/top.depth;
	mid.depth = -1/mid.depth;
	bot.depth = -1/bot.depth;
	split.depth = -1/split.depth;

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
				if (points[c].depth > depths[newX][newY]) {

					depths[newX][newY] = points[c].depth;
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
				if (points[c].depth > depths[newX][newY]) {

					depths[newX][newY] = points[c].depth;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(newX, newY, set);
				}
			}
			

		}

	}

}

void drawCornellWireframe(DrawingWindow &window, std::vector<ModelTriangle> &triangles) {
	for (int i = 0; i < triangles.size(); i++) {
		CanvasTriangle triangle;
		for (int j = 0; j < 3; j++) {
			glm::vec3 cameraToVertex = glm::vec3(triangles[i].vertices[j].x - camera.x, triangles[i].vertices[j].y - camera.y, triangles[i].vertices[j].z - camera.z);

			glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;

			int u = -(distance * (adjustedVector.x)/(adjustedVector.z)) + (window.width / 2);
			int v = (distance * (adjustedVector.y)/(adjustedVector.z)) + (window.height / 2);

			triangle.vertices[j] = CanvasPoint(u, v);
		}
		
		drawTriangle(window, triangle, Colour(255,255,255));
 	}

	glm::vec3 cameraToVertex = glm::vec3(light.x - camera.x, light.y - camera.y, light.z - camera.z);

	glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;

	int u = -(distance * (adjustedVector.x)/(adjustedVector.z)) + (window.width / 2);
	int v = (distance * (adjustedVector.y)/(adjustedVector.z)) + (window.height / 2);

	// prints red pixels to show light location
	window.setPixelColour(u, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	window.setPixelColour(u+1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	window.setPixelColour(u, v+1, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	window.setPixelColour(u-1, v, (255 << 24) + (255 << 16) + (0 << 8) + 0);
	window.setPixelColour(u, v-1, (255 << 24) + (255 << 16) + (0 << 8) + 0);


}

void drawCornell(DrawingWindow &window, std::vector<ModelTriangle> &triangles) {
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

		if (isTexture == true) {
			textureFill(window, triangle, texture, depths);
		} else fillCornell(window, triangle, triangles[i].colour, depths);
		
 	}

}

void raytraceCornell(DrawingWindow &window, std::vector<ModelTriangle> &triangles) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			glm::vec3 falo((WIDTH/2) - x, y - (HEIGHT/2), distance);
			glm::vec3 ray = camera - falo;
			ray = normalize(cameraOrientation * ray);
			RayTriangleIntersection intersect = getClosestIntersection(triangles, ray);
			if (!std::isinf(intersect.distanceFromCamera)) {
				if (inShadow(triangles, intersect.intersectionPoint, intersect.triangleIndex) && shadows) {
					float brightness = 0.11;
					Colour colour = triangles[intersect.triangleIndex].colour;
					colour.red *= brightness;
					colour.blue *= brightness;
					colour.green *= brightness;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(x, y, set);
				} else {
					float brightness = getBrightness(intersect.intersectionPoint, triangles[intersect.triangleIndex].normal);
					Colour colour = triangles[intersect.triangleIndex].colour;
					colour.red *= brightness;
					colour.blue *= brightness;
					colour.green *= brightness;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(x, y, set);
				}

			} 
		}
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
				triangle.normal = glm::cross(glm::vec3(triangle.vertices[2] - triangle.vertices[0]), glm::vec3(triangle.vertices[1] - triangle.vertices[0]));
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
		else if (event.key.keysym.sym == SDLK_b) light.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_g) light.y += 0.1;
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
		else if (event.key.keysym.sym == SDLK_0) {
			resetCamera();
		} 
		else if (event.key.keysym.sym == SDLK_1) {
			drawing = WIREFRAME;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			drawing = RASTERISE;
		}
		else if (event.key.keysym.sym == SDLK_3) {
			drawing = RAYTRACE;
		}
		else if (event.key.keysym.sym == SDLK_LEFTBRACKET) proximity = (proximity) ? false : true;
		else if (event.key.keysym.sym == SDLK_RIGHTBRACKET) angle_of = (angle_of) ? false : true;
		else if (event.key.keysym.sym == SDLK_HASH) shadows = (shadows) ? false : true;
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
		draw(window);
		
		if (drawing == WIREFRAME) {
			drawCornellWireframe(window, triangles);
		} else if (drawing == RASTERISE) {
			drawCornell(window, triangles);
		} else {
			raytraceCornell(window, triangles);
		}

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
