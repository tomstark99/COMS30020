#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240

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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float interval = (to - from) / (numberOfValues-1);
	std::vector<float> result(numberOfValues);
	for(size_t i = 0; i < result.size(); i++) {
		result[i] = from+(interval*i);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<float> intervals(3);
	for(size_t i = 0; i < intervals.size(); i++) {
		intervals[i] = (to[i] - from[i]) / (numberOfValues-1);
	}
	std::vector<glm::vec3> result(numberOfValues);
	for(size_t i = 0; i < result.size(); i++) {
		result[i] = glm::vec3(from[0]+(intervals[0]*i),from[1]+(intervals[1]*i),from[2]+(intervals[2]*i));
	}
	return result;
}

void draw_interpolate(DrawingWindow &window) {
	window.clearPixels();
	std::vector<float> values = interpolateSingleFloats(255, 0, window.width);
	for (size_t y = 0; y < window.width; y++) {
		float greyscale = values[y];
		uint32_t colour = (255 << 24) + (int(greyscale) << 16) + (int(greyscale) << 8) + int(greyscale);
		for (size_t x = 0; x < window.height; x++) {
			window.setPixelColour(y, x, colour);
		}
	}
}

void draw_interpolate_rgb(DrawingWindow &window) {
	window.clearPixels();

	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std::vector<glm::vec3> values_left = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
	std::vector<glm::vec3> values_right = interpolateThreeElementValues(topRight, bottomRight, window.height);

	for (size_t x = 0; x < window.height; x++) {
		std::vector<glm::vec3> values = interpolateThreeElementValues(values_left[x], values_right[x], window.width);
		for (size_t y = 0; y < window.width; y++) {
			uint32_t colour = (255 << 24) + (int(values[y].x) << 16) + (int(values[y].y) << 8) + int(values[y].z);
			window.setPixelColour(y, x, colour);
		}
	}
}

int main(int argc, char *argv[]) {

	// std::vector<glm::vec3> result;
	// glm::vec3 from(1, 4, 9.2);
	// glm::vec3 to(4, 1, 9.8);
	// result = interpolateThreeElementValues(from, to, 4);
	// for(size_t i=0; i<result.size(); i++) std::cout << glm::to_string(result[i]) << " ";
	// std::cout << std::endl;


	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	DrawingWindow window_rgb = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		if (window_rgb.pollForInputEvents(event)) handleEvent(event, window_rgb);
		update(window_grey);
		update(window_rgb);
		draw_interpolate(window_grey);
		draw_interpolate_rgb(window_rgb);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
		window_rgb.renderFrame();
	}
}
