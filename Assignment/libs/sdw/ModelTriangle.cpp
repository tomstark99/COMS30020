#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), texturePoints(), colour(std::move(trigColour)), normal(), mirror() {}

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &n0, const glm::vec3 &n1, const glm::vec3 &n2, Colour trigColour) :
		vertices({{v0, v1, v2}}), normals({{n0, n1, n2}}), texturePoints(), colour(std::move(trigColour)), normal(), mirror() {}

std::ostream &operator<<(std::ostream &os, const ModelTriangle &t) {
	// os << "(" << triangle.normals[0].x << ", " << triangle.normals[0].y << ", " << triangle.normals[0].z << ")\n";
	// os << "(" << triangle.normals[1].x << ", " << triangle.normals[1].y << ", " << triangle.normals[1].z << ")\n";
	// os << "(" << triangle.normals[2].x << ", " << triangle.normals[2].y << ", " << triangle.normals[2].z << ")\n";
	os << "Triangle with colour " << t.colour << " at\n";
	os << "Vertex (1): [" << t.vertices[0].x << "," << t.vertices[0].y << "," << t.vertices[0].z << "]\n";
	os << "Vertex (2): [" << t.vertices[1].x << "," << t.vertices[1].y << "," << t.vertices[1].z << "]\n";
	os << "Vertex (3): [" << t.vertices[2].x << "," << t.vertices[2].y << "," << t.vertices[2].z << "]\n\n";
	os << "Has TexturePoints\n";
	os << "TexturePoint (1): [" << t.texturePoints[0].x << "," << t.texturePoints[0].y << "]\n";
	os << "TexturePoint (2): [" << t.texturePoints[1].x << "," << t.texturePoints[1].y << "]\n";
	os << "TexturePoint (3): [" << t.texturePoints[2].x << "," << t.texturePoints[2].y << "]\n\n";
	return os;
}
