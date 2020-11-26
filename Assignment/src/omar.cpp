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

 //glm::vec3 camera(0.0, 0.0, 4.0);
glm::vec3 camera(140.0,135.0,210.0);
float distance = 700;
glm::mat3 cameraOrientation(
	glm::vec3(1.0, 0.0, 0.0),
	glm::vec3(0.0, 1.0, 0.0),
	glm::vec3(0.0, 0.0, 1.0)
);
int drawing = 1;
glm::vec3 light(0.0,0.85,0.0);
float lightStrength = 30;
bool basic = false;
bool gouraurdDraw = false;
bool phongDraw = true;
//glm::vec3 light(-0.64901096, 2.739334, 0.532032);
//glm::vec3 light(-0.64901096, 2.7384973, -0.51796794);
//glm::vec3 light(0.650989, 2.7384973, -0.51796794);
//glm::vec3 light(0.650989, 2.739334, 0.532032);

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

			if ((u >= 0.0f) && (u <= 1.0f) && (v >= 0.0f) && (v <= 1.0f) && ((u + v) <= 1.0f) && t > 0.05f && t < length) {
				shadow = true;
				break;
			}
		}

	}

	return shadow;
}

float getBrightness(glm::vec3 intersectionPoint, glm::vec3 normal) {
	glm::vec3 lightRay = light - intersectionPoint;
	glm::vec3 cameraRay = (camera * cameraOrientation) - intersectionPoint;
	float length = glm::length(lightRay);

	float angleOfIncidence = glm::dot(glm::normalize(lightRay), normal);
	float brightness = lightStrength/(4 * PI * length*length);

	glm::vec3 angleOfReflection = glm::normalize(lightRay) - ((2.0f*normal)*glm::dot(glm::normalize(lightRay), normal));

	float specular = std::pow(glm::dot(glm::normalize(angleOfReflection), glm::normalize(cameraRay)), 128);

	if (angleOfIncidence > 0) {
		brightness *= angleOfIncidence;
	} 

	if (specular) {
		brightness += specular;
	}

	if (brightness > 1) {
		brightness = 1;
	} 
	if (brightness < 0.2) {
		brightness = 0.2;
	}

	return brightness;
}

float gouraurd(RayTriangleIntersection intersection) {
	glm::vec3 lightRay = light - intersection.intersectionPoint;
	glm::vec3 cameraRay = (camera * cameraOrientation) - intersection.intersectionPoint;
	ModelTriangle triangle = intersection.intersectedTriangle;
	float length = glm::length(lightRay);
	std::vector<float> brightnesses;

	for(int i = 0; i < 3; i++) {
		float temp = glm::dot(triangle.normals[i], glm::normalize(lightRay));
		brightnesses.push_back(temp);
	}
	float angleOfIncidence = (1 - intersection.u - intersection.v) * brightnesses[0] + intersection.u * brightnesses[1] + intersection.v * brightnesses[2];

	float brightness = lightStrength*angleOfIncidence/(4 * PI * length*length);
	// if (angleOfIncidence > 0) {
	// 	brightness *= angleOfIncidence;
	// }

	std::vector<glm::vec3> reflections;
	for(int i = 0; i < 3; i++) {
		glm::vec3 temp = glm::normalize(lightRay) - ((2.0f*triangle.normals[i])*glm::dot(glm::normalize(lightRay), triangle.normals[i]));
		reflections.push_back(temp);
	}
	glm::vec3 angleOfReflection = (1 - intersection.u - intersection.v) * reflections[0] + intersection.u * reflections[1] + intersection.v * reflections[2];

	float specular = std::pow(glm::dot(glm::normalize(angleOfReflection), glm::normalize(cameraRay)), 128);

	if (specular) {
		brightness += specular;
	} 
	if (brightness > 1) {
		brightness = 1;
	} 
	if (brightness < 0.2) {
		brightness = 0.2;
	}

	return brightness;
}

float phong(RayTriangleIntersection intersection) {
	glm::vec3 lightRay = light - intersection.intersectionPoint;
	ModelTriangle triangle = intersection.intersectedTriangle;
	glm::vec3 cameraRay = (camera * cameraOrientation) - intersection.intersectionPoint;
	float length = glm::length(lightRay);
	glm::vec3 specLightRay = intersection.intersectionPoint - light;
	
	glm::vec3 interpolatedNormal = (1 - intersection.u - intersection.v) * triangle.normals[0] + intersection.u * triangle.normals[1] + intersection.v * triangle.normals[2];
	//glm::vec3 interpolatedNormal = triangle.normals[0] + (intersection.u * (triangle.normals[1] - triangle.normals[0])) + (intersection.v * (triangle.normals[2] - triangle.normals[0]));

	float angleOfIncidence = glm::dot(glm::normalize(lightRay), glm::normalize(interpolatedNormal));
	//float brightness = lightStrength*angleOfIncidence/(4 * PI * length*length);
	float brightness = lightStrength/(4 * PI * length*length);

	glm::vec3 angleOfReflection = glm::normalize(specLightRay) - (2.0f*glm::normalize(interpolatedNormal)*glm::dot(glm::normalize(specLightRay), glm::normalize(interpolatedNormal)));

	float specular = std::pow(glm::dot(glm::normalize(angleOfReflection), glm::normalize(cameraRay)), 128);

	if (angleOfIncidence > 0) {
		brightness *= angleOfIncidence;
	} else {
		brightness *= 0;
	}

	if (specular >= 0) {
		brightness += specular*0.2;
	}

	if (brightness > 1) {
		brightness = 1;
	} 
	if (brightness < 0.2) {
		brightness = 0.2;
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
				closestIntersection.u = u;
				closestIntersection.v = v;
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

void texture_half_triangle(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window, std::vector<std::vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];
	
	std::vector<CanvasPoint> left = interpolatePoints(top, mid, abs(mid.y-top.y)+2);
	std::vector<CanvasPoint> right = interpolatePoints(top, bot, abs(mid.y-top.y)+2);
	std::vector<TexturePoint> left_texture = interpolatePoints(top.texturePoint, mid.texturePoint, abs(mid.y-top.y)+2);
	std::vector<TexturePoint> right_texture = interpolatePoints(top.texturePoint, bot.texturePoint, abs(mid.y-top.y)+2);

	for(int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);

		std::vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);
		std::vector<TexturePoint> points_texture = interpolatePoints(left_texture[i], right_texture[i], steps+2);
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

void textureFill(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture, std::vector<std::vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot, mid);
	}
	if (mid.y < top.y) {
		std::swap(mid, top);
	}
	if (bot.y < mid.y) {
		std::swap(bot, mid);
	}

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	split.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	split.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	split.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	CanvasTriangle t_1 = CanvasTriangle(top,mid,split);
	CanvasTriangle t_2 = CanvasTriangle(bot,mid,split);

	texture_half_triangle(t_1, texture, window, depths);
	texture_half_triangle(t_2, texture, window, depths);
}

void fillCornell(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depths) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot, mid);
	}
	if (mid.y < top.y) {
		std::swap(mid, top);
	}
	if (bot.y < mid.y) {
		std::swap(bot, mid);
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

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		CanvasTriangle t;
		for(int j = 0; j < triangle.vertices.size(); j++) {
			glm::vec3 vertex = triangle.vertices[j];
			glm::vec3 cam_to_vertex(vertex.x - camera.x, vertex.y - camera.y, vertex.z - camera.z);
			glm::vec3 adjusted_vertex = cam_to_vertex * cameraOrientation;
		
			int u = -(distance * (adjusted_vertex.x)/(adjusted_vertex.z)) + (window.width/2);
			int v = (distance * (adjusted_vertex.y)/(adjusted_vertex.z)) + (window.height/2);
			;
			t.vertices[j] = CanvasPoint(u,v, adjusted_vertex.z);
			t.vertices[j].texturePoint = triangle.texturePoints[j];
		}
		if (triangle.colour.name != "") {
			TextureMap texture("texture_logo.ppm");
			for(int j = 0; j < t.vertices.size(); j++) {
				t.vertices[j].texturePoint.x *= texture.width;
				t.vertices[j].texturePoint.y *= texture.height;
			}
			std::cout << "ALIEEE" << std::endl;
			textureFill(window, t, texture, depths);
		} else {
			fillCornell(window, t, triangle.colour, depths);
		}
	}
}

	// for (int i = 0; i < triangles.size(); i++) {
	// 	CanvasTriangle triangle;
		// bool isTexture = false;
		// TextureMap texture;

		// if (triangles[i].colour.name != "") {
		// 	texture = TextureMap(triangles[i].colour.name);
		// 	std::cout << triangles[i].colour.name << std::endl;
		// 	isTexture = true;
		// }
		// for (int j = 0; j < 3; j++) {
		// 	glm::vec3 vertex = triangles[i].vertices[j];
		// 	glm::vec3 cameraToVertex = glm::vec3(vertex.x - camera.x, vertex.y - camera.y, vertex.z - camera.z);

		// 	glm::vec3 adjustedVector = cameraToVertex * cameraOrientation;

		// 	int u = -(distance * (adjustedVector.x)/(adjustedVector.z)) + (window.width / 2);
		// 	int v = (distance * (adjustedVector.y)/(adjustedVector.z)) + (window.height / 2);

		// 	triangle.vertices[j] = CanvasPoint(u, v, adjustedVector.z);

		// 	triangle.vertices[j].texturePoint = triangles[i].texturePoints[j];

		// 	// if (isTexture == true) {
		// 	// 	triangle.vertices[j].texturePoint.x *= texture.width;

		// 	// }	
		// }
		// if(triangles[i].colour.name != "") {
		// 	TextureMap texture(triangles[i].colour.name);
		// 	for(int j = 0; j < triangle.vertices.size(); j++) {
		// 		triangle.vertices[j].texturePoint.x *= texture.width;
		// 		triangle.vertices[j].texturePoint.y *= texture.height;
		// 	}
		// 	textureFill(window, triangle, texture, depths);
		// } else fillCornell(window, triangle, triangles[i].colour, depths);

		// if (isTexture == true) {
		// 	std::cout << "myyyyyman" << std::endl;
		
//  	}

// }

void raytraceCornell(DrawingWindow &window, std::vector<ModelTriangle> &triangles) {
	for (int y = 0; y < window.height; y++) {
		for (int x = 0; x < window.width; x++) {
			glm::vec3 falo((WIDTH/2) - x, y - (HEIGHT/2), distance);
			glm::vec3 ray = camera - falo;
			ray = normalize(cameraOrientation * ray);
			RayTriangleIntersection intersect = getClosestIntersection(triangles, ray);
			if (!std::isinf(intersect.distanceFromCamera)) {
				float brightness;
				if (phongDraw) {
					brightness = phong(intersect);
				} else if (basic) {
					brightness = getBrightness(intersect.intersectionPoint, triangles[intersect.triangleIndex].normal);
				} else {
					brightness = gouraurd(intersect);
				}
				bool shadow = inShadow(triangles, intersect.intersectionPoint, intersect.triangleIndex);
				if (shadow) brightness = 0.2;
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

void vertexNormals(std::vector<ModelTriangle> &triangles) {
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle t = triangles[i];
        std::vector<glm::vec3> normals;
        for(int v = 0; v < t.vertices.size(); v++) {
            glm::vec3 vertex = t.normal;
            int count = 1;
            for(int j = 0; j < triangles.size(); j++) {
                ModelTriangle t_ = triangles[j];
                for(int u = 0; u < t_.vertices.size(); u++) {
                    if(i != j && t.vertices[v].x == t_.vertices[u].x && t.vertices[v].y == t_.vertices[u].y && t.vertices[v].z == t_.vertices[u].z) {
                        if (acos(dot(t.normal, t_.normal)/(length(t.normal)*length(t_.normal))) < PI/4) {
							vertex = vertex + t_.normal;
							count = count + 1;
						}
                    }
                }
            }
            vertex = vertex / float(count);
            triangles[i].normals[v] = normalize(vertex);
        }
    }
}

std::vector<ModelTriangle> parseObj(std::string filename, float scale, std::unordered_map<std::string, Colour> colours) {
	std::vector<ModelTriangle> output;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> textureVertices;
	std::string colour;
	std::vector<glm::vec3> normalVecs;

	std::ifstream File(filename);
	std::string line;

	if (filename == "logo.obj") colour = "texture";

	std::cout << colours[colour].name << std::endl;
	
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
				triangle.normal = glm::normalize(glm::cross(glm::vec3(triangle.vertices[1] - triangle.vertices[0]), glm::vec3(triangle.vertices[2] - triangle.vertices[0])));
				if (!normalVecs.empty()) {
					triangle.normals[0] = normalVecs[stoi(a[0])-1];
					triangle.normals[1] = normalVecs[stoi(b[0])-1];
					triangle.normals[2] = normalVecs[stoi(c[0])-1];
				} 
				output.push_back(triangle);
			} else {
				ModelTriangle triangle(vertices[stoi(a[0])-1], vertices[stoi(b[0])-1], vertices[stoi(c[0])-1], colours[colour]);
				//std::cout << triangle.colour << std::endl;
				//std::cout << colours[colour].name << std::endl;
				 //std::cout << colour << std::endl;
				triangle.texturePoints[0] = textureVertices[stoi(a[1])-1];
				triangle.texturePoints[1] = textureVertices[stoi(b[1])-1];
				triangle.texturePoints[2] = textureVertices[stoi(c[1])-1];
				output.push_back(triangle);
			}
		} else if (tokens[0] == "usemtl") {
			colour = tokens[1];
			std::cout << colour << std::endl;
		} else if (tokens[0] == "vt") {
			TexturePoint temp = TexturePoint(stof(tokens[1]), stof(tokens[2]));
			textureVertices.push_back(temp);
		} 
		else if (tokens[0] == "vn") {
			glm::vec3 temp = glm::vec3(stof(tokens[1]), stof(tokens[2]), stof(tokens[3])); 
			normalVecs.push_back(temp);
		}
	}

	if (normalVecs.empty()) {
		vertexNormals(output);
	}

	std::cout << textureVertices.size() << std::endl;

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
			// Colour temp = colours[colour];
			// temp.name = tokens[1];
			// std::cout << temp.name << std::endl;
			// colours[colour] = temp;
			// std::cout << colour << std::endl;
			Colour temp(255, 255, 255);
			temp.name = tokens[1];
			colours.insert({"texture", temp});
			std::cout << colours["texture"] << std::endl;
		}
	}

	File.close();

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
		if (event.key.keysym.sym == SDLK_LEFT) camera.x += 2;
		else if (event.key.keysym.sym == SDLK_RIGHT) camera.x -= 2;
		else if (event.key.keysym.sym == SDLK_UP) camera.y -= 2;
		else if (event.key.keysym.sym == SDLK_DOWN) camera.y += 2;
		else if (event.key.keysym.sym == SDLK_s) camera.z += 0.1;
		else if (event.key.keysym.sym == SDLK_w) camera.z -= 0.1;
		else if (event.key.keysym.sym == SDLK_b) light.y -= 2;
		else if (event.key.keysym.sym == SDLK_n) light.y += 2;
		else if (event.key.keysym.sym == SDLK_z) light.x -= 2;
		else if (event.key.keysym.sym == SDLK_x) light.x += 2; 
		else if (event.key.keysym.sym == SDLK_c) light.z -= 2;
		else if (event.key.keysym.sym == SDLK_v) light.z += 2;
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
		else if (event.key.keysym.sym == SDLK_4) {
			phongDraw = false;
			basic = true;
			gouraurdDraw = false;
			std::cout << "basic" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_5) {
			phongDraw = false;
			basic = false;
			gouraurdDraw = true;
			std::cout << "gouraurd" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_6) {
			phongDraw = true;
			basic = false;
			gouraurdDraw = false;
			std::cout << "phong" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_7) {
			lightStrength -= 1;
			std::cout << "Light Strength: " << lightStrength << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_8) {
			lightStrength += 1;
			std::cout << "Light Strength: " << lightStrength << std::endl;
		}

	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	std::vector<ModelTriangle> triangles;
	std::vector<ModelTriangle> triangles2;
	std::unordered_map<std::string, Colour> colours;

	// colours = parseMtl("cornell-box.mtl");
	// triangles = parseObj("cornell-box.obj", 0.4, colours);

	// triangles2 = parseObj("sphere.obj", 0.4, colours);

	// triangles.insert(triangles.end(), triangles2.begin(), triangles2.end());

	colours = parseMtl("materials.mtl");
	std::cout << colours["texture"] << std::endl;
	triangles = parseObj("logo.obj", 0.3, colours);
	std::cout << triangles.size() << std::endl;

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
