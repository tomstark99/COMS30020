# COMS30020: Computer Graphics Coursework

A Project to demonstrate various modelling, rendering lighting and animation effects on a projection of 3D points.

### Modelling

There are various 3D `.obj` files included with this project, from a sphere and a simple cornell box and hackspace logo to an increasing complex environment combining these models.

`cornell-box.obj` contains the simple cornell box

`cornell-shapes-<num>.obj` contains the cornell box with an increasing number of other models as the appendix number increases

There are also some other misc cornell-box files used during development as well as the low-poly sphere and hackspace logo.

Each `.obj` file uses the same `.mtl` materials file for ease of use.

### Rendering

There are 3 different rendering options when running the program

- wireframe, draws a white outline around all the triangles in the model, creating a mesh. This rendering approach also displays the position of the light within the model
- rasterise, uses a line based interpolation approach to 'fill' triangles giving them a colour.
- raytrace, uses rays projected from the camera to each point in 3D space to find the intersected triangle to display on the 2D plane

### Lighting

- proximity lighting, lights up surfaces with a scale depending on the distance to the light
- angle of incidence lighting, lights up surfaces depending on the angle at which the light hits it (i.e. higher angle means more light spread)
- specular lighting, given a scale factor, makes surfaces appear more shiny or more matte
- ambient lighting, lights up surfaces due to reflections from light off other surfaces

Ambient lighting is computationally expensive therefore a threshold minimum value of `0.3` illuminates the area even when no other lighting is turned on (e.g. from outside sources)

## Usage

### Build and run

To build and run this project use the included `Makefile`
```
make
```
(recommended) or for better results
```
make speedy
```

For animation uncomment the following line inside `main()`

```
animate(t, window_grey, texture);
```

### Keyboard Controls

`1` `2` `3` to switch between raytracing, rasterising and wireframe rendering respectively

`4` `5` `6` to switch between no shading, gourad shading and phong shading respectively

`p` to toggle light position when in raytracing mode (legacy, incompatible with soft shadow approach)

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

`r` to reset the camera view and light position
