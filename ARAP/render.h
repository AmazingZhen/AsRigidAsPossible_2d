#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <unordered_set>

#include "GL/glew.h"
#include "GLFW/glfw3.h"

#include "arap.h"

class Render {
	public:
		void start();

		void selectVertex(float x, float y);
		void moveVertexSelected(float x, float y);

	private:
		bool init();

		void generateSquareMesh();
		void render();

		int findClosestVertexIndex(float x, float y);

	private:
		GLFWwindow* window;

		std::vector<float> originalVertices;
		std::vector<float> vertices;  // {x, y ,z, ..., x, y, z}
		std::vector<unsigned int> indices;

		std::unordered_set<unsigned int> vSelected;

		ARAP solver;
};