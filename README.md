# COMS30020: Computer Graphics

Welcome to the Computer Graphics Unit for 2020

<a href="https://github.com/drslock/CG2020/tree/master/Weekly%20Workbooks">Weekly Workbooks</a> will be released incrementally throughout the teaching block.

## Week 1: RedNoise

Template project to test `SDL`

### usage

run with `make`

## Week 2: Interpolate

Project to demsonstrate pixels and drawing. Output shows two windows

- interpolation of grey values (from 0 to 255) along the x axis and displaying this in an output along the whole y axis
- two dimensional colour interpolation starting with 4 solid colours in 4 corners and interpolating between them, first the left most and right most columns and taking this as the values to interpolate along the x axis

### usage

run with `make`

## Week 3: Triangles

Project to demonstrate line drawing and subsequently triangle drawing (drawing lines between three `CanvasPoint`s), filling and texture mapping.

### usage

run with `make`

`u` draws random colourful outlined triangles within the drawing window

`f` draws random colourful filled triangles within the drawing window

`t` draws a triangle with the texture from `texture.ppm` mapped onto it

## Week 4: Wireframe

Project to demonstrating loading of `.obj` and `.mtl` files and projecting them on a 2D image plane using rasterising

- first with a wireframe using a previous project's `draw_triangle` function
- then by filling the triangles, also using a previous project's `fill_triangle` function.
- then by projecting the triangles with the correct depth by using a depth buffer, taking into account the distance from the camera

### usage

run with `make`

## Week 5: Camera

Project to demonstrate navigation and transformation of a 3D space projected on a 2D plane

- first the floor of the cornell box is textured using the same `texture.ppm` from [week 3](#Week-3:-Triangles)
- being able to move the camera position using translation
- being able to move the camera position using rotation
- orbiting the cornell box while having the camera always look at the origin

### usage

run with `make`

`w` `a` `s` `d` to move forwards, left, back and right respectively within the cornell box

`page_up` `page_down` to move up and down respectively

`q` `e` to rotate around y axis left and right respectively

`0` `9` to rotate around x axis up and down respectively

`left` `right` `up` `down` arrow keys to move camera (look) left, right, up and down respectively

`z` `x` to rotate around z axis (up and down) respectively

`l` to look at origin at any time

`o` to start orbiting (in y axis) around origin

`r` to reset the camera view

## Week 6: Ray

Project to demonstrate the raytracing method of projecting 3D points onto a 2D image plane. For the whole image plane rays are projected out and the closest triangle interaction with that ray is found and the colour is set

Additionally shadows are now possible from setting a light source, again from every point on the 2D image plane, a 3D ray is cast from the surface to the light and any interaction with another triangle means a shadow should be cast

### usage

run with `make` or `make speedy` for better results

`1` `2` `3` to switch between raytracing, rasterising and wireframe rendering respectively

_these are not recommended when in the raytracing rendering method, this is best achieved when in rendering method `2` or `3`_

`w` `a` `s` `d` to move forwards, left, back and right respectively within the cornell box

`page_up` `page_down` to move up and down respectively

`q` `e` to rotate around y axis left and right respectively

`0` `9` to rotate around x axis up and down respectively

`left` `right` `up` `down` arrow keys to move camera (look) left, right, up and down respectively

`z` `x` to rotate around z axis (up and down) respectively

`l` to look at origin at any time

`o` to start orbiting (in y axis) around origin

`r` to reset the camera view

## Week 7: Light

Project to demonstrate various lighting effects on a raytraced projection of 3D points.
- proximity lighting, lights up surfaces with a scale depending on the distance to the light
- angle of incidence lighting, lights up surfaces depending on the angle at which the light hits it (i.e. higher angle means more light spread)
- specular lighting, given a scale factor, makes surfaces appear more shiny or more matte
- ambient lighting, lights up surfaces due to reflections from light off other surfaces

Ambient lighting is computationally expensive therefore a threshold minimum value of `0.3` illuminates the area even when no other lighting is turned on (e.g. from outside sources)

### usage

run with `make` or `make speedy` for better results

`1` `2` `3` to switch between raytracing, rasterising and wireframe rendering respectively

`NUM_8` `NUM_2` `NUM_4` `NUM_6` to move the light forwards, backwards, left and right respectively

`NUM_+` `NUM_-` to move the light up and down respectively

_note. ambient lighting is always on_

`[` to toggle proximity lighting

`]` to toggle angle of incidence lighting

`#` to toggle shadows

`'` to toggle specular lighting

_these are not recommended when in the raytracing rendering method, this is best achieved when in rendering method `2` or `3`_

`w` `a` `s` `d` to move forwards, left, back and right respectively within the cornell box

`page_up` `page_down` to move up and down respectively

`q` `e` to rotate around y axis left and right respectively

`-` `+`(`=`) to rotate around x axis up and down respectively

`left` `right` `up` `down` arrow keys to move camera (look) left, right, up and down respectively

`z` `x` to rotate around z axis (up and down) respectively

`l` to look at origin at any time

`o` to start orbiting (in y axis) around origin

`r` to reset the camera view
