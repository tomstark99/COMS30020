#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
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

void texture_line(CanvasPoint from, CanvasPoint to, TextureMap texture, DrawingWindow &window) {
	float x_diff = to.x - from.x;
	float y_diff = to.y - from.y;
	float x_diff_texture = to.texturePoint.x - from.texturePoint.x;
	float y_diff_texture = to.texturePoint.y - from.texturePoint.y;

	float steps = std::max(abs(x_diff),abs(y_diff));

	float x_step = x_diff/steps;
	float y_step = y_diff/steps;
	float x_step_texture = x_diff_texture/steps;
	float y_step_texture = y_diff_texture/steps;

	for(float i = 0.0; i < steps; i++) {
		float x = from.x + (x_step * i);
		float y = from.y + (y_step * i);
		float x_texture = from.texturePoint.x + (x_step_texture * i);
		float y_texture = from.texturePoint.y + (y_step_texture * i);
		window.setPixelColour(round(x), round(y), texture.pixels[round(y_texture)*texture.width + round(x_texture)]);
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

void fill_half_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];

	float x1_diff = mid.x - top.x;
	float y1_diff = mid.y - top.y;
	float x2_diff = bot.x - top.x;
	float y2_diff = bot.y - top.y;

	float steps1 = std::max(abs(x1_diff),abs(y1_diff));
	float steps2 = std::max(abs(x2_diff),abs(y2_diff));

	float x1_step = x1_diff/steps1;
	float x2_step = x2_diff/steps2;
	float y1_step = y1_diff/steps1;
	float y2_step = y2_diff/steps2;

	for(float i = 0.0; i < steps1; i++) {
		for(float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1_step * i);
			float y1 = top.y + (y1_step * i);
			float x2 = top.x + (x2_step * j);
			float y2 = top.y + (y2_step * j);
			draw_line(CanvasPoint(round(x1),round(y1)),CanvasPoint(round(x2),round(y2)), colour, window);
		}
	}
}

void texture_half_triangle(CanvasTriangle triangle, TextureMap texture, DrawingWindow &window) {
	CanvasPoint top = triangle.vertices[0];
    CanvasPoint mid = triangle.vertices[1];
    CanvasPoint bot = triangle.vertices[2];

	float x1_diff = mid.x - top.x;
	float y1_diff = mid.y - top.y;
	float x2_diff = bot.x - top.x;
	float y2_diff = bot.y - top.y;
	float x1_diff_texture = mid.texturePoint.x - top.texturePoint.x;
	float y1_diff_texture = mid.texturePoint.y - top.texturePoint.y;
	float x2_diff_texture = bot.texturePoint.x - top.texturePoint.x;
	float y2_diff_texture = bot.texturePoint.y - top.texturePoint.y;

	float steps1 = std::max(abs(x1_diff),abs(y1_diff));
	float steps2 = std::max(abs(x2_diff),abs(y2_diff));

	float x1_step = x1_diff/steps1;
	float x2_step = x2_diff/steps2;
	float x1_step_texture = x1_diff_texture/steps1;
	float x2_step_texture = x2_diff_texture/steps2;
	float y1_step = y1_diff/steps1;
	float y2_step = y2_diff/steps2;
	float y1_step_texture = y1_diff_texture/steps1;
	float y2_step_texture = y2_diff_texture/steps2;

	for(float i = 0.0; i < steps1; i++) {
		for(float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1_step * i);
			float y1 = top.y + (y1_step * i);
			float x2 = top.x + (x2_step * j);
			float y2 = top.y + (y2_step * j);

			float texture_x1 = top.texturePoint.x + (x1_step_texture * i);
			float texture_y1 = top.texturePoint.y + (y1_step_texture * i);
			float texture_x2 = top.texturePoint.x + (x2_step_texture * j);
			float texture_y2 = top.texturePoint.y + (y2_step_texture * j);

			CanvasPoint p_1(round(x1),round(y1));
			CanvasPoint p_2(round(x2),round(y2));

			p_1.texturePoint = TexturePoint(round(texture_x1),round(texture_y1));
			p_2.texturePoint = TexturePoint(round(texture_x2),round(texture_y2));
			texture_line(p_1, p_2, texture, window);
		}
	}
}

void fill_triangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
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

	// CanvasPoint mid_2 = find_mid(top, mid, bot);

	// CanvasTriangle t_1 = CanvasTriangle(top,mid,mid_2);
	// CanvasTriangle t_2 = CanvasTriangle(mid_2,mid,bot);

	fill_half_triangle(triangle, colour, window);
	// fill_half_triangle(t_2, colour, window);

	draw_triangle(triangle, Colour(255,255,255), window);
}

void texture_triangle(TextureMap texture, CanvasTriangle triangle, DrawingWindow &window) {
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

	// CanvasPoint mid_2;
	// mid_2.y = mid.y;
	// mid_2.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	// float scale = (mid.y - top.y)/(bot.y-top.y);
	// mid_2.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	// mid_2.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

	// CanvasTriangle t_1 = CanvasTriangle(top,mid,mid_2);
	// CanvasTriangle t_2 = CanvasTriangle(mid_2,mid,bot);

	texture_half_triangle(triangle, texture, window);
	// texture_half_triangle(t_2, texture, window);
	
	draw_triangle(triangle, Colour(255,255,255), window);
}

void random_triangle(DrawingWindow &window, int event) {
	CanvasPoint p_1 = CanvasPoint(std::rand()%window.width,std::rand()%window.height);
	CanvasPoint p_2 = CanvasPoint(std::rand()%window.width,std::rand()%window.height);
	CanvasPoint p_3 = CanvasPoint(std::rand()%window.width,std::rand()%window.height);

	Colour colour = Colour(std::rand()%256, std::rand()%256, std::rand()%256);

	if(event == SDLK_u) {
		draw_triangle(CanvasTriangle(p_1,p_2,p_3), colour, window);
	}
	if(event == SDLK_f) {
		fill_triangle(CanvasTriangle(p_1,p_2,p_3), colour, window);
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u || event.key.keysym.sym == SDLK_f) random_triangle(window, event.key.keysym.sym);
		else if (event.key.keysym.sym == SDLK_t) {
			TextureMap texture("src/texture.ppm");

			CanvasPoint p_1(160,10);
			CanvasPoint p_2(300,230);
			CanvasPoint p_3(10,150);
			p_1.texturePoint = TexturePoint(195,5);
			p_2.texturePoint = TexturePoint(395,380);
			p_3.texturePoint = TexturePoint(65,330);

			CanvasTriangle t(p_1,p_2,p_3);
			texture_triangle(texture, t, window_grey);
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {

	// TextureMap texture("src/texture.ppm");
    // std::cout << texture.width << std::endl;
    // std::cout << texture.height << std::endl;

	// CanvasPoint p_1(160,10);
	// CanvasPoint p_2(300,230);
	// CanvasPoint p_3(10,150);
	// p_1.texturePoint = TexturePoint(195,5);
	// p_2.texturePoint = TexturePoint(395,380);
	// p_3.texturePoint = TexturePoint(65,330);

	// CanvasTriangle t(p_1,p_2,p_3);

	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		update(window_grey);
		// texture_triangle(texture, t, window_grey);
		// fill_triangle(t, Colour(0,0,255),window_grey);
		// draw_triangle(t, Colour(255,255,255), window_grey);
		// random_triangle(window_grey);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}
