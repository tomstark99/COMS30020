#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>

#define WIDTH 320
#define HEIGHT 240

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

void fillHalf(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	CanvasPoint top = triangle[0];
	CanvasPoint mid = triangle[1];
	CanvasPoint bot = triangle[2];

	float x1Diff = mid.x - top.x;
	float y1Diff = mid.y - top.y;

	float x2Diff = bot.x - top.x;
	float y2Diff = bot.y - top.y;

	float steps1 = std::max(abs(x1Diff), abs(y1Diff));
	float steps2 = std::max(abs(x2Diff), abs(y2Diff));

	float x1Step = x1Diff/steps1;
	float x2Step = x2Diff/steps2;

	float y1Step = y1Diff/steps1;
	float y2Step = y2Diff/steps2;

	for (float i = 0.0; i < steps1; i++) {
		for (float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1Step * i);
			float y1 = top.y + (y1Step * 1);

			float x2 = top.x + (x2Step * j);
			float y2 = top.y + (y2Step * j);

			CanvasPoint point1(round(x1), round(y1));
			CanvasPoint point2(round(x2), round(y2));
			
			drawLine(window, point1, point2, colour);
		}
	}
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	// getting flat coords

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	CanvasTriangle triangle1(top, mid, split);
	CanvasTriangle triangle2(split, mid, bot);

	std::cout << top << std::endl;
	std::cout << mid << std::endl;
	std::cout << split << std::endl;

	fillHalf(window, triangle1, colour);
	fillHalf(window, triangle2, colour);
	drawTriangle(window, triangle, Colour(255,255,255));

}

void draw(DrawingWindow &window) {
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
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl; 
		else if (event.key.keysym.sym == SDLK_u) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand()%window.width, rand()%window.height), 
				CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));

			Colour colour = Colour(rand()%256, rand()%256, rand()%256);
			drawTriangle(window, triangle, colour); 
		}
		else if (event.key.keysym.sym == SDLK_f) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand()%window.width, rand()%window.height), 
				CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));

			Colour colour = Colour(rand()%256, rand()%256, rand()%256);
			drawFilledTriangle(window, triangle, colour); 
		}
		
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}

// Sort into vertical order
// Divide triangle into 2
