#include "ModelTriangle.h"
#include <utility>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), texturePoints(), colour(std::move(trigColour)), normal() {}

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &n0, const glm::vec3 &n1, const glm::vec3 &n2, Colour trigColour) :
		vertices({{v0, v1, v2}}), normals({{n0, n1, n2}}), texturePoints(), colour(std::move(trigColour)), normal() {}

std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.normals[0].x << ", " << triangle.normals[0].y << ", " << triangle.normals[0].z << ")\n";
	os << "(" << triangle.normals[1].x << ", " << triangle.normals[1].y << ", " << triangle.normals[1].z << ")\n";
	os << "(" << triangle.normals[2].x << ", " << triangle.normals[2].y << ", " << triangle.normals[2].z << ")\n";
	return os;
}
